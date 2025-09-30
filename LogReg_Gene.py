# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 14:47:42 2025

@author: Nayara
"""

#%% Pacotes 

import pandas as pd
import numpy as np 
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor

#%% Upload de dados transcriptomicos

data_gene = pd.read_csv('readCounts_placenta_updated.csv')
data_gene.info()
data_gene.set_index('GeneID', inplace= True)
data_gene.head()

metadados_gene = pd.read_csv('metadados_placenta_updated.csv')
metadados_gene.info()
metadados_gene.rename(columns ={'Unnamed: 0' : 'GSM'}, inplace = True)
metadados_gene.drop('source_name', axis = 1, inplace = True)
metadados_gene.head()

### Genes de interesse
genes_degs = pd.read_excel('Genes_final.xlsx')
genes_degs.head()

map_genes = genes_degs['GeneID']
filt_genes_counts = data_gene.loc[data_gene.index.isin(map_genes)]

#%% Normalização e escalonamento dos genes

### Normalização LogCPM
cpm = filt_genes_counts/(filt_genes_counts.sum(axis=0))*1e6
log_cpm = np.log2(cpm + 1).T


### Escalonamento - Z Scale
z_logcpm = StandardScaler().fit_transform(log_cpm)

df_zlogcpm = pd.DataFrame(z_logcpm,
                          columns= log_cpm.columns,
                          index = log_cpm.index)

### Separação de somente os genes de interesse - DEGs
dic_genes = dict(zip(genes_degs['GeneID'], genes_degs['Gene'] ))
df_zlogcpm.rename(columns = dic_genes, inplace = True)
dic_group = dict(zip(metadados_gene['GSM'],metadados_gene['Group']))
df_zlogcpm['group'] = df_zlogcpm.index.map(dic_group)

#%% Upload de dados clínicos

clinical = pd.read_csv('Dados_IMC_cat.csv')
clinical.head()

clinical_data = clinical.loc[:, ['GSM', 'group', 'ogtt_1h', 
       'maternal_weight_adjusted', 'BMI_adjust',
        'group_dummie', 'BMI_adjust_cat', 'Sobrepeso']]


#Padronizaçao dos nomes das colunas
clinical_data.rename(columns={'ogtt_1h' :'oGTT_1h', 'maternal_weight_adjusted':'Peso',
                            'BMI_adjust':'IMC',
                            'group_dummie': 'Grupos'}, inplace = True)

#%% Trabalhando com dados Escalonados

### Genes de interesse 

data = df_zlogcpm.copy()
data.reset_index(inplace = True)
data.rename(columns = {'index' :'GSM'}, inplace =True)

genes_interesse = data.loc[:, ['GSM', 'CSH1',
                               'CSH2','CSHL1',
                               'GH2','GH1', 'group']]

### Merge dos dados Clinicos com os genes de interesse
all_data = pd.merge(clinical_data,
                    genes_interesse[['GSM', 'CSH1','CSH2','CSHL1','GH2', 'GH1']],
                    on='GSM',how = 'left')
all_data.info()

#%% Imputação MICE para os dados escalonados 

# Dicionário para armazenar os dataframes imputados de cada grupo
imputed_dfs = {}

# Iterar sobre cada grupo para realizar a imputação separadamente
for group in all_data['group'].unique():
    
    # Filtrar o dataframe para o grupo atual
    group_df = all_data[all_data['group'] == group].copy()
    
    # Selecionar apenas as colunas numéricas relevantes para a imputação
    imputation_cols = ['oGTT_1h','Peso','IMC',
                       'CSH1','CSH2','CSHL1','GH2', 'GH1']
    imputation_data = group_df[imputation_cols]
    
    # Criar uma instância do imputador MICE
    imputer = IterativeImputer(random_state=42)
    
    # Imputar os valores
    imputed_data = imputer.fit_transform(imputation_data)
    
    # Converter a matriz de volta para um dataframe
    imputed_df = pd.DataFrame(imputed_data, columns=imputation_cols, index=group_df.index)
    
    # Acessar os valores originais para obter o mínimo e máximo 
    original_csh2 = all_data.loc[all_data['group'] == group, 'CSH2'].dropna()
    original_csh1= all_data.loc[all_data['group'] == group, 'CSH1'].dropna()
    original_cshl1= all_data.loc[all_data['group'] == group, 'CSHL1'].dropna()
    original_gh2= all_data.loc[all_data['group'] == group, 'GH2'].dropna()
    original_gh1= all_data.loc[all_data['group'] == group, 'GH1'].dropna()
   
    # Ajustar os valores imputados para respeitar o mínimo e máximo
    imputed_df['CSH2'] = imputed_df['CSH2'].clip(lower=original_csh2.min(), upper=original_csh2.max())
    imputed_df['CSH1'] = imputed_df['CSH1'].clip(lower=original_csh1.min(), upper=original_csh1.max())
    imputed_df['CSHL1'] = imputed_df['CSHL1'].clip(lower=original_cshl1.min(), upper=original_cshl1.max())
    imputed_df['GH2'] = imputed_df['GH2'].clip(lower=original_gh2.min(), upper=original_gh2.max())
    imputed_df['GH1'] = imputed_df['GH1'].clip(lower=original_gh1.min(), upper=original_gh1.max())
  
   
    # Substituir as colunas originais pelas colunas imputadas
    group_df.loc[:, imputation_cols] = imputed_df
    
    # Armazenar o dataframe imputado
    imputed_dfs[group] = group_df
    
# Concatenar os dataframes de volta em um único dataframe final
df_imputed = pd.concat(imputed_dfs.values()).sort_index()

# Testando se os dados foram imputados
df_imputed.loc[:, ['CSH1','CSH2','CSHL1','GH2', 'GH1']].isna().value_counts()
df_imputed.loc[:, ['CSH1','CSH2','CSHL1','GH2', 'GH1']].describe()

#%% HEATMAP 

# Calcular as médias por grupo
df_means = df_imputed.groupby('group')[['CSH1','CSH2','GH2']].mean()

# Mapear os nomes dos grupos para português
df_means.index = df_means.index.map({'Control': 'Control', 'GDM': 'GDM'})

# Criar heatmap com as médias dos grupos
plt.figure(figsize=(2.5, 2.5), dpi=300) 
heatmap = sns.heatmap(df_means,
                     cmap='coolwarm',  
                     center=0,
                     square=True,     
                     cbar_kws={'label': '', 'shrink': 0.6, 'aspect': 15},  
                     linewidths=0.3,   
                     linecolor='white',
                     annot=False,      
                     cbar=True)        

# Remove labels e título
plt.ylabel('')
plt.xlabel('')
plt.title('')

# Formata os ticks - tamanhos reduzidos para figura pequena
plt.yticks(fontsize=8, fontweight='bold', rotation=0)
plt.xticks(fontsize=8, fontweight='bold', rotation=45)

# Formata a barra de cores
cbar = heatmap.collections[0].colorbar
cbar.ax.tick_params(labelsize=6)  # Ticks menores para barra
cbar.set_label('Expression mean (Z-score)', size=6)        # Label vazia para barra

# Ajustar layout de forma muito compacta
plt.tight_layout(pad=0.5)  # Padding mínimo

plt.show()


#%%  Score de genes transcriptomicos 

df_imputed['CSH_group'] = df_imputed[['CSH2','CSH1']].sum(axis=1)

#%% Modelos Utilizando dados Escalonados 

### MODELO 1 - oGTT, Peso e Genes  
modelo1 = sm.Logit.from_formula('Grupos ~ oGTT_1h + Peso + CSH_group',
                     data=df_imputed).fit()
#AIC / BIC
aic_m1= modelo1.aic
bic_m1 = modelo1.bic
print(modelo1.summary())
print(f"Modelo 1 : AIC({aic_m1:.2f}) e BIC({bic_m1:.2f})")

print (100*'-')

### MODELO 2 - oGTT, IMC e Genes 
modelo2 = sm.Logit.from_formula('Grupos ~ oGTT_1h + IMC + CSH_group',
                    data=df_imputed).fit()
#AIC / BIC
aic_m2= modelo2.aic
bic_m2 = modelo2.bic
print(modelo2.summary())
print(f"Model 2 : AIC({aic_m2:.2f}) e BIC({bic_m2:.2f})")

print (100*'-')


### MODELO 3 - oGTT
modelo3 = sm.Logit.from_formula('Grupos ~ oGTT_1h + Sobrepeso + CSH_group',
                    data=df_imputed).fit()
#AIC / BIC
aic_m3= modelo3.aic
bic_m3 = modelo3.bic
print(modelo3.summary())
print(f"Model 3 : AIC({aic_m3:.2f}) e BIC({bic_m3:.2f})")

#%% Colinearidade entre variaveis com Genes Escalonados 

### Matriz de correlação 
correlation_matrix = df_imputed[['oGTT_1h',
                                 'Peso',
                                 'IMC',
                                 'CSH_group']].corr()

print(correlation_matrix)

plt.figure(figsize=(8, 6))
sns.heatmap(correlation_matrix,
            annot=True,
            cmap='coolwarm',
            fmt=".2f", 
            linewidths=.5)
plt.show()



##VIF - Variance Inflation Factor

X_VIF = df_imputed[['oGTT_1h','Peso','CSH_group']]
X_VIF_const = sm.add_constant(X_VIF)

VIF = pd.DataFrame()
VIF["variable"] = X_VIF_const.columns
VIF["VIF"] = [variance_inflation_factor(X_VIF_const.values, i)
              for i in range(X_VIF_const.shape[1])]
VIF["Tolerance"] = 1 / VIF["VIF"]

print(VIF)


X_VIF = df_imputed[['oGTT_1h','IMC','CSH_group']]
X_VIF_const = sm.add_constant(X_VIF)

VIF = pd.DataFrame()
VIF["variable"] = X_VIF_const.columns
VIF["VIF"] = [variance_inflation_factor(X_VIF_const.values, i)
              for i in range(X_VIF_const.shape[1])]
VIF["Tolerance"] = 1 / VIF["VIF"]

print(VIF)

#%% ODDs RATIOs 

def calcular_or(model):
    params = model.params
    conf = model.conf_int()
    conf.columns = ['2.5%', '97.5%']

    or_values = np.exp(params)
    or_conf = np.exp(conf)
    
    results = pd.DataFrame({
        'Coeficiente': params,
        'Razão de Chaces': or_values,
        'IC 2.5%': or_conf['2.5%'],
        'IC 97.5%': or_conf['97.5%'],
        'p-value': model.pvalues})
    
    return results

or_m1 = calcular_or(modelo1).round(3)
or_m2 = calcular_or(modelo2).round(3)
or_m3 = calcular_or(modelo3).round(3)

print("\n=== Modelo 1 ===")
print(or_m1.T)
print("\n=== Modelo 2 ===")
print(or_m2.T)
print("\n=== Modelo 3 ===")
print(or_m3.T)


#%% Exportando tabela dados imputados 

#df_imputed.to_excel('DadosClinicosGenes.xlsx', index = False)

#%% Teste de interação IMC - CSH

#### Centralizar a variável IMC na média
df_imputed['IMC_c'] = df_imputed['IMC'] - df_imputed['IMC'].mean()
df_imputed['CSH_group_c'] = df_imputed['CSH_group'] - df_imputed['CSH_group'].mean()
df_imputed['oGTT_1h_c'] = df_imputed['oGTT_1h'] - df_imputed['oGTT_1h'].mean()


modelo4 = sm.Logit.from_formula('Grupos ~ oGTT_1h_c + IMC_c + CSH_group_c + CSH_group_c:oGTT_1h_c',
                                     data=df_imputed).fit()

print(modelo4 .summary())


# Seleciona variáveis e renomeia
correlation_matrix = df_imputed[['oGTT_1h', 'IMC', 'CSH_group']].rename(
    columns={'oGTT_1h': 'Glicemia-1h', 'CSH_group' : 'CSH'}
).corr()

print(correlation_matrix)

# Plot
plt.figure(figsize=(2, 2), dpi=300)
sns.heatmap(correlation_matrix,
            annot=True,
            cmap='coolwarm',
            fmt=".2f",
            linewidths=.5)

# Ajusta rotação dos labels
plt.xticks(rotation=45, ha="right", fontsize=8)
plt.yticks(rotation=360, ha="right", fontsize=8)

plt.show()

