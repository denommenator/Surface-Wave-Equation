#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 00:02:25 2020

@author: robertdenomme

Generate the adjacency matrix of the meshgrid.

Theory is that for NxN grid with basis indexed by pairs (r,s) 0<=r<n 
in lexicographical order, the pattern is a_ij=1 iff |i-j|=1, N. The 1
comes from the first element adjacency, and the N from the second element
adjacency
"""
import numpy as np
import scipy.linalg
from tqdm import tqdm

X = np.array([[1,2,3],[4,5,6],[7,8,9],[10,11,12]], dtype=float)
X_dot = np.zeros(np.shape(X))

#practice data with zero velocity
    
class WaveEngine:
    """This class defines functions for working with the wave equation on
    a rectangular grid. the displacements are organized in an MxN matrix
    X, and the finite difference X(t)-X(t-1) by X_dot. This package defines
    phase space as a single vector, given by indexing the coordinates of
    X and X_dot, then concatenating the two. The culmination is the definition
    of a dynamics matrix D of shape (2M*N, 2M*N) which acts on a phases space
    vector."""
    
    def __init__(self, shape, Nfrm, alpha=.01, fps=.1, delta_x = 1):
        #numerical fields
        self.shape = shape
        self.alpha=alpha
        self.fps = fps
        self.delta_t = fps
        self.delta_x = delta_x
        
        #matrix fields that are calculated upon init
        self.L = self.get_laplacian_matrix()
        self.D_cont = self.get_dynamics_matrix()
        
        
        self.Nfrm = Nfrm
    
    def get_phase_position_index(self,i,j):
        """More of an internal helper function to help transiton to phase space"""
        M, N = self.shape
        return i*N+j
    
    def get_config_index(self,p):
        """More of an internal helper function to help transiton to phase space"""
        M, N = self.shape
        return (p//N, p%N)
    
    def flatten(self,X):
        """More of an internal helper function to help transiton to phase space"""
        M, N = self.shape
        V = np.ndarray((M*N), dtype=float)
        for i,dum in enumerate(X):
            for j,dum in enumerate(X[0]):
                V[self.get_phase_position_index(i,j)]=X[i][j]
        return V
    
    def unflatten(self,V):
        """internale helper function if useful"""
        M, N = self.shape
        X = np.ndarray((M,N),dtype=float)
        for p,value in enumerate(V):
            (i,j)=self.get_config_index(p)
            X[i][j]=value
        return X
    
    def to_phase_space(self,X,X_dot):
        """Send the pair X,X_dot of matrices to a single vector in order to use
        matrix notation for the dynamic map X(t)->X(t+1)"""
        M, N = self.shape
        return np.concatenate((self.flatten(X),self.flatten(X_dot)),0)
    
    def from_phase_space(self,P):
        """Bring a single phase space vector back into the pair (X,X_dot) of 
        matrices. Note: all dynamics happen in phase space"""
        M, N = self.shape
        X=np.ndarray((M,N),dtype=float)
        X_dot=np.ndarray((M,N),dtype=float)
        
        
        for p,dum in enumerate(P):
            if(p<M*N):
                (i,j)=self.get_config_index(p)
                X[i][j]=P[p]
            else:
                (i,j)=self.get_config_index(p-M*N)
                X_dot[i][j]=P[p]
        return (X,X_dot)
    
    
    
    def is_adjacent(self,p,p_prime):
        (i,j)=self.get_config_index(p)
        (i_prime,j_prime)=self.get_config_index(p_prime)
        if i==i_prime and j-j_prime in {1,-1}:
            return True
        elif i-i_prime in {1,-1} and j==j_prime:
            return True
        else:
            return False
    
    def get_laplacian_matrix(self):
        """Make the M*NxM*N lagrangian matrix corresponding to the (first half 
        of the) phase space 
        coordinates. The form is A-4I, where A adds the 4 neighbors of i,j 
        together"""
        M, N = self.shape
        L = np.ndarray((M*N,M*N),dtype=float)
        for p in range(M*N):
            for p_prime in range(M*N):
                (i,j) = self.get_config_index(p)
                (i_prime,j_prime)=self.get_config_index(p_prime)
                if p==p_prime:
                    L[p,p_prime]=-4
                elif self.is_adjacent(p,p_prime) is True:
                    L[p,p_prime]=1
                else:
                    L[p,p_prime]=0
        return L
                            
    def get_discrete_dynamics_matrix(self):
        """This is it. This computes the dynamics of the system 
        X(t+1) = X(t)+(X(t)-X(t-1))+alpha*L*X(t)
        in phase space. Note X_dot(t) = X(t)-X(t-1), which is the second
        coordinate of phase space. From here you can get the spectral 
        qualities of the system"""
        M, N = self.shape
        alpha = self.alpha
        L=self.get_laplacian_matrix()
        D = np.ndarray((2*M*N,2*M*N),dtype=float)
        I = np.identity(M*N)
        # Z = np.zeros((M*N,M*N))
        #These values were for the older version ith (X(t),X(t+1)) instead of (X)
        # row_one = np.concatenate((Z,I),axis=1)
        # row_two = np.concatenate((-1*I,2*I+alpha*self.get_lagrangian_matrix()),axis=1)
        
        row_one = np.concatenate((I+alpha*L,I),axis=1)
        row_two = np.concatenate((alpha*L,I),axis=1)
        
        
        D = np.concatenate((row_one,row_two),axis=0)
        return D
    
    def get_dynamics_matrix(self):
        """This returns a matrixc D_cont for which 
        <X(i*delta_t),X'(i*delta_t)> = D_cont^i*<X(0),X'(0)>"""
        M, N = self.shape
        alpha = self.alpha
        delta_t=self.delta_t
        L=self.L
        
        
        D_log = np.ndarray((2*M*N,2*M*N),dtype=float)
        I = np.identity(M*N, dtype=float)
        Zero = 0*I
        
        row_one = np.concatenate((Zero,I),axis=1)
        row_two = np.concatenate((alpha**2*L,Zero),axis=1)
        D_log = np.concatenate((row_one,row_two),axis=0)
        D_cont = scipy.linalg.expm(delta_t*D_log)
        
        return D_cont
    
    def get_continuous_trajectory(self, X, X_dot):
        M, N = self.shape
        Nfrm = self.Nfrm
        P = self.to_phase_space(X,X_dot)
        D_cont = self.D_cont
        
        Z_phase = [P]
        for i in range(1,Nfrm):
            """sorry for the weird indexing. Want to multiply by D_cont 
            every step except the first one"""
            Zi_phase = D_cont@Z_phase[-1]
            Z_phase.append(Zi_phase)
        
        Z = [self.from_phase_space(Zi_phase)[0] for Zi_phase in Z_phase]
        return Z
    
    def get_next_from_phase(self,P):
        """This should be re-implemented using convolutions since 
        that is way faster"""
        return self.D_cont@P
    
    
    def get_next(self,X,X_dot):
        """Find and return the next pair from the dynamics of the 
        wave equation. This is the main method for animating the wave"""
        P = self.to_phase_space(X,X_dot)
        P_next = self.get_next_from_phase(P)
        return self.from_phase_space(P_next)
        
        

