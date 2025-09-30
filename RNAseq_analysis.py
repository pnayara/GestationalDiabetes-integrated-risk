# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 11:24:09 2025

@author: Nayara
"""

#%% Pacotes
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns 
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms


#%%  Reading tsv 

NCBI_file = pd.read_csv('GSE203346_raw_counts_GRCh38.p13_NCBI.tsv', sep='\t')
print(NCBI_file.head())
print(NCBI_file.info())
print(NCBI_file.columns)
print(NCBI_file.shape)
   
    
#%% Metadados 
metadata_df = pd.read_csv('metadata.csv')
print(metadata_df.head())
print(metadata_df.info())
print(metadata_df.columns)

    
#%% Verificar a compatibilidade entre os dados do NCBI e Metadados


NCBI_samples = NCBI_file.columns
NCBI_samples = NCBI_samples.drop(NCBI_samples[0])
equal_test = set(NCBI_samples) == set(metadata_df['Sample Name'])
print('Todos os valores estão presentes', equal_test)

out_in_NCBI = set(metadata_df['Sample Name']) - set(NCBI_samples) 
print("Amostras em metadata_df que não estão em NCBI_samples:", out_in_NCBI)

samples_out = ['GSM6167415', 'GSM6167412', 'GSM6167423', 'GSM6167356']
samples_out_meta = metadata_df[metadata_df['Sample Name'].isin(samples_out)]

#%% Alinhamento Metadados - NCBI

# Verificar a distribuição dos grupos e tecidos 
print("Contagem de amostras por 'Group':")
print(metadata_df['Group'].value_counts())
print("\nContagem de amostras por 'source_name':")
print(metadata_df['source_name'].value_counts())

# Selecionar colunas de interesse e definir 'Sample Name' como índice
colunas_interesse = ['Sample Name', 'Group', 'source_name']
metadata_subset_df = metadata_df[colunas_interesse].copy() 

# Definir 'Sample Name' como o índice do DataFrame de metadados
metadata_subset_df.set_index('Sample Name', inplace=True)
print(metadata_subset_df.head())
print(metadata_subset_df.shape)

#Alinhar os metadados à matriz de contagens
NCBI_file= NCBI_file.set_index('GeneID')
sample_order_in_counts = NCBI_file.columns

# Reordenar o metadata_subset_df para seguir essa ordem
# .loc[sample_order_in_counts] seleciona as linhas na ordem especificada
metadata_aligned_df = metadata_subset_df.loc[sample_order_in_counts]
print(metadata_aligned_df.head())
print(metadata_aligned_df.shape)

# Verificação final: Assegurar que os índices dos metadados alinhados
# correspondem exatamente às colunas da matriz de contagens
if list(metadata_aligned_df.index) == list(NCBI_file.columns):
    print("\nOs nomes e a ordem das amostras nos metadados alinhados correspondem às colunas da matriz de contagens!")
else:
    print("\nA ordem/nomes das amostras não correspondem. Investigar!")

#avaliando grupos     
metadata_aligned_df.groupby('Group')['source_name'].value_counts()

#%%  Dataframes 

NCBI_file # read_table
metadata_aligned_df # metadados alinhados 

print(NCBI_file.shape)
print(metadata_aligned_df.shape)


#%% Separar os SETS de Amostras

#Placenta
metadados_placenta = metadata_aligned_df[metadata_aligned_df['source_name'
                                                             ] == 'Placenta']
counts_placenta = NCBI_file[metadados_placenta.index]
print(counts_placenta.shape)

counts_placenta  #read_table
metadados_placenta # metadados placenta


#Umbilical cord Blood
metadados_umbilical = metadata_aligned_df[metadata_aligned_df['source_name'
                                                             ] == 'Umbilical cord blood']
counts_umbilical = NCBI_file[metadados_umbilical.index]
print(counts_umbilical .shape)

counts_umbilical #read_table
metadados_umbilical # Metadados Umbilical 

#%% Library Size

lib_size_placenta = counts_placenta.sum(axis = 0).to_frame(name='LibrarySize')
lib_size_placenta = lib_size_placenta.join(metadados_placenta[['Group']], how ='inner')
print(lib_size_placenta.head())
print(lib_size_placenta.groupby('Group').describe().map(lambda x: f"{x:0.1f}").T)

lib_size_umbilical = counts_umbilical.sum(axis = 0).to_frame(name='LibrarySize')
lib_size_umbilical = lib_size_umbilical.join(metadados_umbilical[['Group']], how ='inner')
print(lib_size_umbilical.head())


#%% PLOTs Library Size

plt.figure(figsize=(16, 8))

plt.subplot(1, 2, 1)
sns.barplot(x= lib_size_placenta.index,
            y='LibrarySize',
            data=lib_size_placenta.sort_values(by=['Group','LibrarySize']),
            hue='Group',
            width = 0.8,
            linewidth =0.5,
            edgecolor = 'black',
            palette=['#ff7f0e','#1f77b4'],
            dodge=False)
plt.title('Library Size - Placenta', fontsize = 18,
          fontweight = 'bold',
          fontname= 'Arial', pad =20)
plt.xlabel('Samples',fontsize = 16,
          fontweight = 'bold',
          fontname= 'Arial',
          labelpad=10)
plt.ylabel('Total Counts (10^7)',
           fontsize = 16,
          fontweight = 'bold',
          fontname= 'Arial',
          labelpad=10)
plt.legend(title='',
           frameon=True,
           title_fontsize=12,
           fontsize=14,
           handlelength=1.0,
           handleheight=1.0 )
plt.xticks(rotation=90,
           ha='right',
           va = 'top',
           fontsize=14,
           fontweight = 'bold',
           fontname= 'Arial',
           rotation_mode = 'anchor')
plt.yticks(fontsize=16,
           fontweight = 'bold',
           fontname= 'Arial')
sns.despine(right=True)
plt.tight_layout()

plt.subplot(1, 2, 2)
sns.barplot(x= lib_size_umbilical.index,
            y='LibrarySize',
            data=lib_size_umbilical.sort_values(by=['Group','LibrarySize']),
            hue='Group',
            width = 0.8,
            linewidth =0.5,
            edgecolor = 'black',
            palette=['#ff7f0e','#1f77b4'],
            dodge=False)
plt.title('Library Size - Umbilical Blood', fontsize = 18,
          fontweight = 'bold',
          fontname= 'Arial', pad =20)
plt.xlabel('Samples',fontsize = 16,
          fontweight = 'bold',
          fontname= 'Arial',
          labelpad=10)
plt.ylabel('Total Counts (10^7)',
           fontsize = 16,
          fontweight = 'bold',
          fontname= 'Arial',
          labelpad=10)
plt.legend(title='',
           frameon=True,
           title_fontsize=12,
           fontsize=14,
           handlelength=1.0,
           handleheight=1.0 )
plt.xticks(rotation=90,
           ha='right',
           va = 'top',
           fontsize=14,
           fontweight = 'bold',
           fontname= 'Arial',
           rotation_mode = 'anchor')
plt.yticks(fontsize=16,
           fontweight = 'bold',
           fontname= 'Arial')
sns.despine(right=True)
plt.tight_layout()

plt.show()

# BOXPLOT

plt.figure(figsize=(12, 10))

plt.subplot(1, 2, 1)
sns.boxplot(x='Group',
            y='LibrarySize', 
            data=lib_size_placenta,
            order = ['Control', 'GDM'],
            hue='Group',
            width = 0.5,
            linewidth = 1,
            boxprops={'linewidth': 1.5,
                      'edgecolor': 'black'},  
            whiskerprops={'linewidth': 1.5}, 
            capprops={'linewidth': 1.5},
            palette=['#1f77b4','#ff7f0e'],
            showmeans = True)
plt.title('Library Size Distribution - Placenta',
          fontsize = 18,
          fontweight = 'bold',
          fontname= 'Arial',
          pad =20)          
plt.ylabel('Total Counts (10^7)',
           fontsize = 16,
           fontweight = 'bold',
           fontname= 'Arial',
           labelpad=10)
plt.xlabel('')
plt.xticks(fontsize=14,
           fontweight = 'bold',
           fontname= 'Arial')
plt.yticks(fontsize=16,
           fontweight = 'bold',
           fontname= 'Arial')
plt.tight_layout()
sns.despine(right=True)


plt.subplot(1, 2, 2)
sns.boxplot(x='Group',
            y='LibrarySize', 
            data=lib_size_umbilical,
            order = ['Control', 'GDM'],
            hue='Group',
            width = 0.5,
            linewidth = 1,
            boxprops={'linewidth': 1.5,
                      'edgecolor': 'black'},  
            whiskerprops={'linewidth': 1.5}, 
            capprops={'linewidth': 1.5},
            palette=['#1f77b4','#ff7f0e'],
            showmeans = True)
plt.title('Library Size Distribution - Umbilical',
          fontsize = 18,
          fontweight = 'bold',
          fontname= 'Arial',
          pad =20)          
plt.ylabel('Total Counts (10^7)',
           fontsize = 16,
           fontweight = 'bold',
           fontname= 'Arial',
           labelpad=10)
plt.xlabel('')
plt.xticks(fontsize=14,
           fontweight = 'bold',
           fontname= 'Arial')
plt.yticks(fontsize=16,
           fontweight = 'bold',
           fontname= 'Arial')
plt.tight_layout()
sns.despine(right=True)


plt.show()

#%% Expressão total por amostra

#Placenta
cpm_placenta = (counts_placenta.div(lib_size_placenta['LibrarySize']))*1e6
log_cpm_placenta = np.log2(cpm_placenta + 1)
expression_placenta = log_cpm_placenta.sum(axis=0).to_frame(name='Total Expression')
expression_placenta = expression_placenta.join(metadados_placenta[['Group']], how = 'inner')

#Umbilical Cord
cpm_umbilical = (counts_umbilical.div(lib_size_umbilical['LibrarySize']))*1e6
log_cpm_umbilical = np.log2(cpm_umbilical + 1)
expression_umbilical = log_cpm_umbilical.sum(axis=0).to_frame(name='Total Expression')
expression_umbilical = expression_umbilical.join(metadados_umbilical[['Group']], how = 'inner')

#%% PLOT Expressão total por amostra

plt.figure(figsize=(12, 10))

plt.subplot(1, 2, 1)
sns.boxplot(x='Group',
            y='Total Expression', 
            data=expression_placenta,
            order = ['Control', 'GDM'],
            hue='Group',
            width = 0.5,
            linewidth = 1,
            boxprops={'linewidth': 1.5,
                      'edgecolor': 'black'},  
            whiskerprops={'linewidth': 1.5}, 
            capprops={'linewidth': 1.5},
            palette=['#1f77b4','#ff7f0e'],
            showmeans = True)
plt.title('Total Expression - Placenta',
          fontsize = 18,
          fontweight = 'bold',
          fontname= 'Arial',
          pad =20)          
plt.ylabel('log2(CPM + 1)',
           fontsize = 16,
           fontweight = 'bold',
           fontname= 'Arial',
           labelpad=10)
plt.xlabel('')
plt.xticks(fontsize=14,
           fontweight = 'bold',
           fontname= 'Arial')
plt.yticks(fontsize=16,
           fontweight = 'bold',
           fontname= 'Arial')
plt.tight_layout()
sns.despine(right=True)


plt.subplot(1, 2, 2)
sns.boxplot(x='Group',
            y='Total Expression', 
            data=expression_umbilical,
            order = ['Control', 'GDM'],
            hue='Group',
            width = 0.5,
            linewidth = 1,
            boxprops={'linewidth': 1.5,
                      'edgecolor': 'black'},  
            whiskerprops={'linewidth': 1.5}, 
            capprops={'linewidth': 1.5},
            palette=['#1f77b4','#ff7f0e'],
            showmeans = True)
plt.title('Total Expression- Umbilical',
          fontsize = 18,
          fontweight = 'bold',
          fontname= 'Arial',
          pad =20)          
plt.ylabel('log2(CPM + 1)',
           fontsize = 16,
           fontweight = 'bold',
           fontname= 'Arial',
           labelpad=10)
plt.xlabel('')
plt.xticks(fontsize=14,
           fontweight = 'bold',
           fontname= 'Arial')
plt.yticks(fontsize=16,
           fontweight = 'bold',
           fontname= 'Arial')
plt.tight_layout()
sns.despine(right=True)


plt.show()


#%% PCA Sem Filtragem dos Genes de baixa expressão 

# Escalonameto Placenta
scal_logcpm_plac = StandardScaler().fit_transform(log_cpm_placenta.T)
#PCA Placenta
pca_placenta = PCA(n_components=2)
pca_placenta_result = pca_placenta.fit_transform(scal_logcpm_plac)
#df PCA Placenta
pca_placenta_df = pd.DataFrame( data = pca_placenta_result,
                               columns =['PC1', 'PC2'],
                               index = log_cpm_placenta.T.index)
pca_placenta_df = pca_placenta_df.join(metadados_placenta[['Group']],
                                       how = 'inner')

# Escalonamento Umbilical 
scal_logcpm_umbilical = StandardScaler().fit_transform(log_cpm_umbilical.T)
#PCA Umbilical
pca_umbilical = PCA(n_components=2)
pca_umbilical_result = pca_umbilical.fit_transform(scal_logcpm_umbilical)
#df PCA Umbilical
pca_umbilical_df = pd.DataFrame( data = pca_umbilical_result,
                               columns =['PC1', 'PC2'],
                               index = log_cpm_umbilical.T.index)
pca_umbilical_df = pca_umbilical_df.join(metadados_umbilical[['Group']],
                                       how = 'inner')

#%% PLOT PCAs  Sem Filtragem dos Genes de baixa expressão 

plt.figure(figsize=(16, 8))
plt.subplot(1, 2, 1)
sns.scatterplot(
        x="PC1",
        y="PC2",
        hue="Group",
        hue_order=['Control', 'GDM'],
        s=80,
        data=pca_placenta_df)

plt.xticks(fontsize =16)
plt.yticks(fontsize =16)
plt.title('PCA Placenta',
          fontsize =16, fontweight = 'bold')
plt.xlabel(f'PC1: {pca_placenta.explained_variance_ratio_[0]*100:.2f}%',
           fontsize =16, fontweight = 'bold')
plt.ylabel(f'PC2: {pca_placenta.explained_variance_ratio_[1]*100:.2f}%', 
           fontsize =16, fontweight = 'bold')
plt.legend(fontsize =14)
plt.grid(False)
sns.despine()

plt.subplot(1, 2, 2)

sns.scatterplot(
        x="PC1",
        y="PC2",
        hue="Group",
        hue_order=['Control', 'GDM'],
        s=80,
        data=pca_umbilical_df)

plt.xticks(fontsize =16)
plt.yticks(fontsize =16)
plt.title('PCA Umbilical',
          fontsize =16, fontweight = 'bold')
plt.xlabel(f'PC1: {pca_umbilical.explained_variance_ratio_[0]*100:.2f}%',
           fontsize =16, fontweight = 'bold')
plt.ylabel(f'PC2: {pca_umbilical.explained_variance_ratio_[1]*100:.2f}%', 
           fontsize =16, fontweight = 'bold')
plt.legend(fontsize =14)
plt.grid(False)
sns.despine()
 
plt.show()

#%% Remoção de Genes com baixa expressão
genes_mantidos_placenta = (cpm_placenta > 1).sum(axis=1) >= 18
cpm_placenta_filtrado = cpm_placenta[genes_mantidos_placenta]
print(cpm_placenta.shape)
print(cpm_placenta_filtrado.shape)


genes_mantidos_umbilical = (cpm_umbilical > 1).sum(axis=1) >= 20
cpm_umbilical_filtrado = cpm_umbilical[genes_mantidos_umbilical]
print(cpm_umbilical.shape)
print(cpm_umbilical_filtrado.shape)

#%% PCA pós filtragem 

#Log
log_cpm_plac_filt = np.log2(cpm_placenta_filtrado + 1)

# Escalonameto Placenta
scal_logcpm_plac_filt = StandardScaler().fit_transform(log_cpm_plac_filt.T)
#PCA Placenta
pca_placenta_filt = PCA(n_components=2)
pca_placenta_result_filt = pca_placenta_filt.fit_transform(scal_logcpm_plac_filt)
#df PCA Placenta
pca_placenta_df_filt = pd.DataFrame(data = pca_placenta_result_filt,
                               columns =['PC1', 'PC2'],
                               index = log_cpm_plac_filt.T.index)
pca_placenta_df_filt = pca_placenta_df_filt.join(metadados_placenta[['Group']],
                                       how = 'inner')
#Log
log_cpm_umbilical_filt = np.log2(cpm_umbilical_filtrado + 1)

# Escalonamento Umbilical 
scal_logcpm_umbilical_filt = StandardScaler().fit_transform(log_cpm_umbilical_filt.T)
#PCA Umbilical
pca_umbilical_filt = PCA(n_components=2)
pca_umbilical_result_filt = pca_umbilical_filt.fit_transform(scal_logcpm_umbilical_filt)
#df PCA Umbilical
pca_umbilical_df_filt = pd.DataFrame( data = pca_umbilical_result_filt,
                               columns =['PC1', 'PC2'],
                               index = log_cpm_umbilical_filt.T.index)
pca_umbilical_df_filt = pca_umbilical_df_filt.join(metadados_umbilical[['Group']],
                                       how = 'inner')
#%% PLOT PCAs Com filtragem 

plt.figure(figsize=(16, 8), dpi =300)
plt.subplot(1, 2, 1)
sns.scatterplot(data=pca_placenta_df_filt,
        x="PC1",
        y="PC2",
        hue="Group",
        hue_order=['Control', 'GDM'],
        s=80,
        palette= ['#ff7f0e','#1f77b4'])

plt.xticks(fontsize =16)
plt.yticks(fontsize =16)
plt.title('PCA Placenta',
          fontsize =16, fontweight = 'bold')
plt.xlabel(f'PC1: {pca_placenta_filt.explained_variance_ratio_[0]*100:.2f}%',
           fontsize =16, fontweight = 'bold')
plt.ylabel(f'PC2: {pca_placenta_filt.explained_variance_ratio_[1]*100:.2f}%', 
           fontsize =16, fontweight = 'bold')
plt.legend(fontsize =14)
plt.grid(False)
sns.despine()

plt.subplot(1, 2, 2)

sns.scatterplot(data=pca_umbilical_df_filt,
        x="PC1",
        y="PC2",
        hue="Group",
        hue_order=['Control', 'GDM'],
        s=80,
        palette=['#ff7f0e','#1f77b4'])
        
plt.xticks(fontsize =16)
plt.yticks(fontsize =16)
plt.title('PCA Umbilical',
          fontsize =16, fontweight = 'bold')
plt.xlabel(f'PC1: {pca_umbilical_filt.explained_variance_ratio_[0]*100:.2f}%',
           fontsize =16, fontweight = 'bold')
plt.ylabel(f'PC2: {pca_umbilical_filt.explained_variance_ratio_[1]*100:.2f}%', 
           fontsize =16, fontweight = 'bold')
plt.legend(fontsize =16)
plt.grid(False)
sns.despine()
 
plt.show()

#%% Estudo dos Outliers

#Distancia Euclidiana ao Centro (0,0)
pca_placenta_result_filt
pca_umbilical_result_filt

dist_plac = np.sqrt(pca_placenta_result_filt[:,0]**2 + 
                    pca_placenta_result_filt[:,1]**2)

dist_umbilical = np.sqrt(pca_umbilical_result_filt[:,0]**2 + 
                         pca_umbilical_result_filt[:,1]**2)

#Criterio de Tukey  para Outliers
q1,q3 = np.percentile(dist_plac,[25,75])
iqr =q3-q1
limite_sup =q3 +1.5*iqr
outliers_plac = dist_plac >limite_sup
outliers_samples_plac = log_cpm_plac_filt.T.index[outliers_plac]
print('Amostras Outliers na Placenta:', outliers_samples_plac )
#Amostras Outliers na Placenta: Index(['GSM6167371', 'GSM6167388']

q1u , q3u = np.percentile(dist_umbilical,[25,75])
iqr_u =q3u - q1u
limite_sup_u =q3u +1.5*iqr_u
outliers_umbilical = dist_umbilical >limite_sup_u
outliers_samples_umbilical = log_cpm_umbilical_filt.T.index[outliers_umbilical]
print('Amostras Outliers no Umbilical:', outliers_samples_umbilical)

#PLOT Placenta - unico com outliers
plt.figure(figsize=(8, 6))
plt.scatter(pca_placenta_result_filt[:, 0],
            pca_placenta_result_filt[:, 1],
            c='grey',
            label='Samples')
plt.scatter(pca_placenta_result_filt[outliers_plac, 0],
            pca_placenta_result_filt[outliers_plac, 1],
            c='red',
            label='Outliers')
plt.xlabel(f'PC1: {pca_placenta_filt.explained_variance_ratio_[0]*100:.2f}%')
plt.ylabel(f'PC2: {pca_placenta_filt.explained_variance_ratio_[1]*100:.2f}%')
plt.title("Outliers via PCA - Placenta")
plt.legend()
plt.show()

#%% Remoção de 2 outliers Placenta 

#Antes de remover as amostras 'GSM6167371', 'GSM6167388'
cpm_placenta_filtrado 
metadados_placenta

metadados_plac_rem = metadados_placenta.drop(['GSM6167371',
                                              'GSM6167388',
                                              'GSM6167424'], axis = 0)

cpm_plac_fil_rem = cpm_placenta_filtrado.drop(['GSM6167371',
                                              'GSM6167388',
                                              'GSM6167424'], axis = 1)
print(cpm_placenta_filtrado.shape)
print(cpm_plac_fil_rem.shape)

print(metadados_placenta.shape)
print(metadados_plac_rem.shape)

#PCA sem os outliers

#Log
logcpm_plac_rem = np.log2(cpm_plac_fil_rem + 1)
# Escalonameto Placenta
scal_logcpm_plac_rem = StandardScaler().fit_transform(logcpm_plac_rem.T)
#PCA Placenta
pca_placenta_rem = PCA(n_components=2)
pca_placenta_result_rem = pca_placenta_rem.fit_transform(scal_logcpm_plac_rem)
#df PCA Placenta
pca_placenta_df_rem = pd.DataFrame(data = pca_placenta_result_rem,
                               columns =['PC1', 'PC2'],
                               index = logcpm_plac_rem.T.index)
pca_placenta_df_rem = pca_placenta_df_rem.join(metadados_plac_rem[['Group']],
                                       how = 'inner')
print(pca_placenta_df_rem['Group'].value_counts())


#PLOT
plt.figure(figsize=(8, 6))
sns.scatterplot(data=pca_placenta_df_rem,
        x="PC1",
        y="PC2",
        hue="Group",
        hue_order=['Control', 'GDM'],
        s=80,
        palette= ['#ff7f0e','#1f77b4'])

plt.xticks(fontsize =16)
plt.yticks(fontsize =16)
plt.title('PCA Placenta -Outlier Removed',
          fontsize =16, fontweight = 'bold')
plt.xlabel(f'PC1: {pca_placenta_rem.explained_variance_ratio_[0]*100:.2f}%',
           fontsize =16, fontweight = 'bold')
plt.ylabel(f'PC2: {pca_placenta_rem.explained_variance_ratio_[1]*100:.2f}%', 
           fontsize =16, fontweight = 'bold')
plt.legend(fontsize =14)
plt.grid(False)
sns.despine()
plt.show()

#Variancias
explained_variance_ratio = pca_placenta_rem.explained_variance_ratio_
cumulative_explained_variance = np.cumsum(explained_variance_ratio)
print(explained_variance_ratio)
print(cumulative_explained_variance)

#%% Elipses de IC

# Função para plotar elipses de confiança (adaptada da documentação do matplotlib)
def confidence_ellipse(x, y, ax, n_std=1.96, **kwargs):
    if x.size != y.size:
        raise ValueError("x e y devem ter o mesmo tamanho")
    
    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, **kwargs)
    
    mean_x, mean_y = np.mean(x), np.mean(y)
    scale_x = np.sqrt(cov[0, 0]) * n_std
    scale_y = np.sqrt(cov[1, 1]) * n_std
    
    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)
    
    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)

# Plot PCA com elipses


palette_order = {"Control": "lightblue", "GDM": "salmon"}
plt.figure(figsize=(10, 7), dpi = 300)
ax = plt.gca()
sns.scatterplot(
        x="PC1",
        y="PC2",
        hue="Group",
        hue_order=['Control', 'GDM'],
        palette = palette_order,
        s=200,
        data=pca_placenta_df_rem,
        ax=ax)

color = ["lightblue" ,  "salmon"]
groups = ['Control', 'GDM']
for i, group in enumerate(groups):
    group_data = pca_placenta_df_rem [pca_placenta_df_rem ["Group"] == group]
    confidence_ellipse(
        x=group_data["PC1"].values,
        y=group_data["PC2"].values,
        ax=ax,
        n_std=1.96,
        alpha=0.2,
        color=color[i])

plt.xticks(fontsize =16)
plt.yticks(fontsize =16)
plt.title('PCA Placenta',
          fontsize =20, fontweight = 'bold', pad = 20)
plt.xlabel(f'PC1: {pca_placenta_rem.explained_variance_ratio_[0]*100:.2f}%',
           fontsize =18, fontweight = 'bold')
plt.ylabel(f'PC2: {pca_placenta_rem.explained_variance_ratio_[1]*100:.2f}%', 
           fontsize =18, fontweight = 'bold')

legend = ax.legend_
if legend is not None:
    legend.get_title().set_text("")
    for text, new_label in zip(legend.get_texts(), ['Control', 'GDM']):
        text.set_text(new_label)
    # Ajusta o tamanho da fonte
    plt.setp(legend.get_texts(), fontsize=16)

plt.grid(False)
sns.despine()
plt.show()

#%% Datasets finais

#cout_tables 
reads_plac_final = counts_placenta.loc[genes_mantidos_placenta].drop(['GSM6167371',
                                              'GSM6167388',
                                              'GSM6167424'], axis = 1)
print('Before Filter:', counts_placenta.shape)
print('After Filter:', reads_plac_final.shape)

#metadados
meta_plac_final = metadados_placenta.drop(['GSM6167371',
                                              'GSM6167388',
                                              'GSM6167424'], axis = 0)

print('Before Filter:', metadados_placenta.shape)
print('After Filter:', meta_plac_final.shape)


#cout_tables Umbilical
reads_umbilical_final = counts_umbilical.loc[genes_mantidos_umbilical]
print('Before Filter:', counts_umbilical.shape)
print('After Filter:', reads_umbilical_final.shape)
# Metadados continua o mesmo já que não houve remoção de amostras
metadados_umbilical
print(metadados_umbilical.shape)

# %% Dataframes Finais após filtragem 
#placenta
meta_plac_final
reads_plac_final

#meta_plac_final.to_csv('metadados_placenta_updated.csv', index = True)
#reads_plac_final.to_csv('readCounts_placenta_updated.csv', index = True)

#umbilical
metadados_umbilical
reads_umbilical_final

#metadados_umbilical.to_csv('metadados_umbilical_updated.csv', index = True)
#reads_umbilical_final.to_csv('readCounts_umbilical_updated.csv', index = True)


