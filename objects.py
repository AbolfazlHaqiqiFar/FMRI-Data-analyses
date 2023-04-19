import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import networkx as nx
import plotly.graph_objs as go


#------------------------------------------------------------------------------------

class Statics:
    def __init__(self, filename) -> None:
        self.filename = filename
        self.data = self.read_data()
        self.correlation = self.calculate_correlation()

    #--------------------------------------------------------------------------------
    def read_data(self):
        with open(self.filename, 'r') as f:
        # Assuming the data is comma-separated (CSV format), but you can adjust this for your file
            data = pd.read_csv(f, delimiter=' ', header=None)
        return data
    
    #------------------------------------------------------------------------------
    def calculate_correlation(self):
        correlation = self.data.corr()
        return correlation
    
    #--------------------------------------------------------------------------------
    def signed_network(self, threshold=0):
        correlation = self.calculate_correlation()
        signnet = correlation.applymap(lambda x: 1 if x > threshold else (-1 if x < -threshold else 0))
        return signnet
    
    #----------------------------------------------------------------------------------------
    def network_hamiltonian(self):
        signnet = self.signed_network()
        n = signnet.shape[0]
        count =0
        network_hamiltonian=0
        for i in range(n-1):
            for j in range(i+1, n):
                count += 1
                network_hamiltonian += signnet.iloc[i,j]
        return network_hamiltonian / count
    
    #-------------------------------------------------------------------------------------------
    def interaction_hamiltonian(self):
        signnet = self.signed_network()
        n = signnet.shape[0]
        interaction_hamiltonian = 0
        count = 0
        for i in range(n-1):
            for j in range(i+1, n):
                for k in range(n-1):
                    for l in range(k+1, n):
                        if i!= k and l!= j:
                            count += 1
                            interaction_hamiltonian += signnet.iloc[i,j]*signnet.iloc[k,l] /2
        return interaction_hamiltonian / count
    
    #---------------------------------------------------------------------------------------------
    def balance_hamiltonian(self):
        signnet = self.signed_network()
        n = signnet.shape[0]
        balance_hamiltonian = 0
        count = 0
        for i in range(n-2):
            for j in range(i+1, n-1):
                for k in range(j+1, n):
                    count +=1
                    balance_hamiltonian += signnet.iloc[i,j]*signnet.iloc[j,k]*signnet.iloc[k,i] 
        return balance_hamiltonian/ count
    
    #----------------------------------------------------------------------------------------------
    def carpet_plot(self, data):
        colorscale = [[0, '#EAECEE'], [1, '#1E90FF']]
        trace = go.Heatmap(z=data.values, colorscale=colorscale)
        layout = go.Layout(title='Carpet plot of thresholded correlation matrix', 
                           xaxis=dict(title='Genes'), yaxis=dict(title='Genes'))
        fig = go.Figure(data=[trace], layout=layout)
        fig.show()
    
    

        