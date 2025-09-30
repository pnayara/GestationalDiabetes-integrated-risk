# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 19:34:56 2025

@author: Nayara
"""

#%% Instalação do pacote que fara integração R - Python

#pip install rpy2
#pip install gseapy


#%% Testando a integração

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter


base = importr('base')
print("Biblioteca 'base' carregada!")

edgeR= importr("edgeR")

limma = importr("limma")
print("Biblioteca 'limma' carregada!")

r_version = base.R_Version()
print(f"Versão do R instalada: {r_version.rx2('version.string')[0]}")


#%% Pacotes Python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp
from gseapy import prerank

#%% Datasets utilzados 

#placenta 
tableCounts_placenta = pd.read_csv('readCounts_placenta_updated.csv')
tableCounts_placenta = tableCounts_placenta.set_index('GeneID')
metadados_placenta = pd.read_csv('metadados_placenta_updated.csv')
metadados_placenta = metadados_placenta.rename(columns ={'Unnamed: 0' : 'GSM'})

print(tableCounts_placenta.head())
print(metadados_placenta.head())

#%% Ativar a conversão automática de pandas <-> R data.frames

with localconverter(ro.default_converter + pandas2ri.converter):
    r_counts_placenta = ro.conversion.py2rpy(tableCounts_placenta)
    r_metadata_placenta = ro.conversion.py2rpy(metadados_placenta)
    
#%% Passar objetos para o ambiente R

ro.globalenv['counts'] = r_counts_placenta
ro.globalenv['metadata'] = r_metadata_placenta

#%% Pipeline limma-voom para análise de DEGs (GDM vs Control) - Integração R

#criar Objeto DGEList
ro.r('dge <- DGEList(counts = counts)')

#normalizar pelos fatores da biblioteca
ro.r('dge <- calcNormFactors(dge)')

#Criar o design matrix
ro.r('''
     group <- factor(metadata$Group, levels =c("Control", 'GDM'))
     design <- model.matrix(~group)
     colnames(design) <- c("Intercept", "GDMvsControl")
     ''')

#aplica voom e ajusta o modelo
ro.r('v <- voom( dge, design, plot = TRUE)')

#ajusta o modelo linear
ro.r('fit <- lmFit (v, design)')

#Moderação empirica de variancia
ro.r('fit <- eBayes(fit)')

#Obter os DEGs
ro.r('deg_results <- topTable(fit, coef ="GDMvsControl", number =Inf, sort.by ="P")')

#%% Trazer o 'deg_results' para o python para trabalhar no pandinha <3

with localconverter(ro.default_converter + pandas2ri.converter):
    deg_results_df = ro.conversion.rpy2py(ro.r['deg_results'])

print(deg_results_df.head())
print(deg_results_df.describe().map(lambda x: f"{x:0.3f}").T)

#%% anotação dos genes - Placenta

annot = pd.read_csv('Human.GRCh38.p13.annot.tsv', sep = '\t')
annot['GeneID'] = annot['GeneID'].astype('str')
print(annot.head())
print(annot.columns)

deg_results_df =deg_results_df.reset_index()
deg_results_df =deg_results_df.rename(columns ={'index' :'GeneID'})

deg_annot_placenta = deg_results_df.merge(annot[['GeneID',
                                        'Symbol',
                                        'Description',
                                        'GeneType']], on= 'GeneID', how= 'left')

#deg_annot_placenta.to_csv('DEGs_placenta_limmaVoon.csv', index = False)


#%% Placenta - Parte 2

mod_group = pd.read_csv('DEGs_placenta_limmaVoon.csv')
print(mod_group.head())



#%% Explorando os Genes DEGs da placenta 

#VULCANO PLOT - No filter

log10= -np.log10(mod_group['P.Value'])
red = (mod_group['logFC']> 1) & (mod_group['P.Value'] < 0.05)
blue = (mod_group['logFC'] < -1) & (mod_group['P.Value'] < 0.05)

colors = np.where(red, 'red', 
                np.where(blue, 'blue', 'gray')) ## Essa estrutura faz loopings condicionais em arrays grande.
                                                    # mais eficiente para muitos dados


plt.figure(figsize=(8, 6))
plt.scatter(mod_group['logFC'], 
            log10, 
            c=colors,
            s=mod_group['AveExpr']*7,  # usa average expression como tamanho
            alpha=0.5,
            edgecolors='k')

plt.axhline(-np.log10(0.05),
            color='gray', 
            linestyle='--',
            dashes = (5,2),
            linewidth=0.5)
plt.axvline(1, 
            color='gray',
            linestyle='-',
            linewidth=0.5)
plt.axvline(-1, 
            color='gray',
            linestyle='-',
            linewidth=0.5)

plt.title('DMG vs Controle',
          fontsize = 18,
          fontweight = 'bold',
          fontname= 'Arial',
          pad =20)
plt.xlabel('Fold change$\mathbf{(log_{2})}$ ',
           fontsize = 16,
           fontweight = 'bold',
           fontname= 'Arial',
           labelpad=10)
plt.ylabel('$\mathbf{-log_{10}}$(p-value)',
           fontsize = 16,
           fontweight = 'bold',
           fontname= 'Arial',
           labelpad=10)
plt.xticks(fontsize = 14, 
           fontname= 'Arial')
plt.yticks(fontsize = 14, 
           fontname= 'Arial')
plt.ylim(0, 5)

plt.grid(False)
sns.despine()
plt.tight_layout()

plt.show()

#%% Analise de vias - preparação

top_genes_up = mod_group[(mod_group['logFC']> 1) & (mod_group['P.Value'] < 0.05)]
top_genes_down = mod_group[(mod_group['logFC'] <-1) & (mod_group['P.Value'] < 0.05)]

#Genes com |LogFC > 1| (up and down), P-value<0.05, AveExpr > 1 e B > -1 - muito conservador
top_genes_filt = mod_group[(abs(mod_group['logFC'])> 1)
                           & (mod_group['P.Value'] < 0.05)
                           & (mod_group['AveExpr'] > 1)
                           & (mod_group['B'] > -1)]

#Filtro um pouco menos conservador: LogFC > 1 (up), P-value<0.05, AveExpr > 1
top_genes_2 = mod_group[(mod_group['logFC']> 1)
                           & (mod_group['P.Value'] < 0.05)
                           & (mod_group['AveExpr'] > 1)]




#%% Estudo dos Genes DEG upregulados

up_genes = top_genes_2['Symbol'].to_list()

# Rodar enrichr para mapeamento contra GO
enrichr_results = gp.enrichr(
    gene_list=up_genes,
    gene_sets='GO_Biological_Process_2021',
    organism='Human',
    outdir=None ) # não salva no disco

# Ver os top 10 processos biológicos
df_enrich = enrichr_results.results
df_enrich[['Term', 'Adjusted P-value', 'Genes']].head(10)


#%% Estudo GSEA 

all_degs = mod_group[['logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B','Symbol',]]
all_degs_filt = all_degs[['logFC', 'Symbol']].sort_values( by = 'logFC', ascending = False)
all_degs_filt = all_degs_filt.set_index('Symbol')['logFC'].sort_values(ascending =False)



prerank_results = prerank(rnk=all_degs_filt,
                          gene_sets='GO_Biological_Process_2021',  
                          outdir=None,
                          permutation_num=1000, # aumentar para maior robustez
                          seed=42)
prerank_results_df = prerank_results.res2d



KEGG_results = prerank(rnk=all_degs_filt,
                          gene_sets='KEGG_2021_Human',  
                          outdir=None,
                          permutation_num=1000, # aumentar para maior robustez
                          seed=42)
KEGG_results_df = KEGG_results.res2d

#%% Criar uma integração entre os genes mais presentes em todas as listas

# Um DataFrame final com os genes da lista exploratória
#Indica presença nas vias (ORA e/ou GSEA)
#Conta quantas vias cada gene aparece
#Atribui uma classificação inicial de prioridade


df_integrativo = top_genes_2[['Symbol','P.Value' ,'logFC','AveExpr']].copy()
df_integrativo.rename(columns ={'Symbol' : 'Gene'}, inplace = True)

#Lista de genes no ORA
genes_in_ora = set()
for gene_list in df_enrich['Genes']:
    genes = gene_list.split(';')
    genes_in_ora.update([g.strip() for g in genes])

#Lista de Genes em GSEA
genes_in_gsea = set()
for gene_list in KEGG_results_df['Lead_genes']:
    genes = gene_list.split(';')
    genes_in_gsea.update([g.strip() for g in genes])

#Verificando genes em comum
df_integrativo['in_ORA'] = df_integrativo['Gene'].isin(genes_in_ora)
df_integrativo['in_GSEA'] = df_integrativo['Gene'].isin(genes_in_gsea)

#Contagem de numero de vias que cada gene aparece
from collections import Counter

ora_gene_counts = Counter()
for gene_list in df_enrich['Genes']:
    genes = [g.strip() for g in gene_list.split(';')]
    ora_gene_counts.update(genes)

gsea_gene_counts = Counter()
for gene_list in KEGG_results_df['Lead_genes']:
    genes = [g.strip() for g in gene_list.split(';')]
    gsea_gene_counts.update(genes)

#Soma de ocorrencia nas duas vias
df_integrativo['N_vias_total'] = df_integrativo['Gene'].apply(
    lambda g : ora_gene_counts.get( g, 0) + gsea_gene_counts.get(g, 0))

#Priorização
def classificar(row):
    in_both = row['in_ORA'] and row['in_GSEA']
    in_either = row['in_ORA'] or row['in_GSEA']
    significant = row['P.Value'] < 0.05
    high_fc = row['logFC'] >= 2
    many_pathways = row['N_vias_total'] >= 5

    # Prioridade 1: Muito Alta 
    if in_both and high_fc and significant and many_pathways:
        return 'Muito Alta'
    if in_either and high_fc and significant and many_pathways:
        return 'Muito Alta'
    
    # Prioridade 2: Alta 
    if in_both and significant and many_pathways:
        return 'Alta'
    if in_either and significant and many_pathways:
        return 'Alta'
    
    # Prioridade 3: 
    if in_both and many_pathways:
        return 'Média'
    if in_either and many_pathways:
        return 'Média'
    
    # Padrão: Baixa
    return 'Baixa'

df_integrativo['avaliacao'] = df_integrativo.apply(classificar, axis = 1)

df_integrativo['avaliacao'].value_counts().plot(kind='bar', title='Distribuição da Prioridade dos Genes')

## Separar os genes Muito alto e alto
genes_selecionados = df_integrativo[df_integrativo['avaliacao'].isin(['Alta', 'Muito Alta'])]
gene_hub = genes_selecionados['Gene']
#gene_hub.to_csv('gene_hub.csv')

#%% Analise de REDE (HUB- Genes)

import networkx as nx

ppi_df = pd.read_csv('string_interactions.tsv', sep ='\t')
ppi_df.info()


#ppi_df é mais completo 
ppi_df.rename(columns={'#node1' : 'node1'}, inplace = True)

#Construindo a Rede 

#Grafo de Interações
grafo = nx.from_pandas_edgelist(ppi_df, source ='node1', target ='node2')
print(f"Num.Nos:{grafo.number_of_nodes()}")
print(f"Num.arestas:{grafo.number_of_edges()}")

#Métricas
degree_dict = dict(grafo.degree())
centrality = nx.betweenness_centrality(grafo)
closeness = nx.closeness_centrality(grafo)

hub_metrics = pd.DataFrame({
    'Gene':list(degree_dict.keys()),
    'Degree': list(degree_dict.values()),
    'Betweenness':[centrality[g] for g in degree_dict.keys()],
    'Closeness' : [closeness[g] for g in degree_dict.keys()]
    })
hub_metrics.sort_values(by='Degree', ascending = False, inplace= True)
hub_metrics.head(10)

## Top 10 Genes com mais conexões
top_hubs = hub_metrics.head(10)
top_hubs = top_hubs['Gene'].tolist()

grafo_2 = grafo.subgraph(top_hubs)

#PLOT
plt.figure(figsize=(8, 6))
pos = nx.spring_layout(grafo_2, seed=42)
nx.draw(grafo_2,
        pos,
        with_labels=True,
        node_size=1200,
        node_color='lightblue', edge_color='gray')
plt.title("Sub-rede dos Top Hub Genes")
plt.show()

#%% Tabela final de Genes 

top_hubs # Genes Centrais 
mod_group # DEG's com todas as informações 
hub_metrics # Dados de Networking
genes_selecionados # Quantidade de vias que aparecem  

df1 = mod_group[mod_group['Symbol'].isin(top_hubs)]
df1.rename( columns ={'Symbol':'Gene'}, inplace=True)

df2 = hub_metrics[hub_metrics['Gene'].isin(top_hubs)]
df3 = genes_selecionados[genes_selecionados['Gene'].isin(top_hubs)]

genes_final =pd.merge(df1, df2 , on='Gene', how='inner')
genes_final = pd.merge(genes_final, df3[['N_vias_total','Gene']], on ='Gene', how = 'inner')

#genes_final.to_csv('Genes_final.csv', index = False )
#genes_final.to_excel('Genes_final.xlsx', index = False )
