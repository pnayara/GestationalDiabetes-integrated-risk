# -*- coding: utf-8 -*-
"""
Created on Mon Sep  8 16:41:24 2025

@author: Nayara
"""


#%% Pacotes 

import pandas as pd
import numpy as np
from scipy.stats import shapiro
from scipy.stats import levene
from scipy.stats import fisher_exact  
import matplotlib.pyplot as plt
import seaborn as sns
import pingouin as pg
import statsmodels.api as sm

#%% Dados

### Dados para os modelos utilizados
# OGTT1h, Peso, IMC, Categorias e Grupo CSH
data = pd.read_excel('DadosClinicosGenes.xlsx') 
data.info()

### Demais dados Clinicos 
# oGTT2h, Idade, Altura
data_2 = pd.read_excel('df_placenta.xlsx', sheet_name= 'data3')
data_2.info()

#%% Preparação dos dados 

### oGTT1h, Peso, IMC, GSH e Categorias
data_stats = data.copy()
control = data_stats.loc[data_stats['group'] == 'Control' ]
gdm = data_stats.loc[data_stats['group'] == 'GDM' ]

### Altura, Idade e oGTT2h
data_stats2 = data_2.copy()
control2 = data_stats2.loc[data_stats2['group'] == 'Control' ]
gdm2 = data_stats2.loc[data_stats2['group'] == 'GDM' ]


#%% Testes Estatísticos 

#### OGTT1h 

#NORMALIDADE
stat_shapiro, pvalue_shapiro = shapiro(control['oGTT_1h'])
print(f"Stat_shapiro:{stat_shapiro:.2f}, p_value:{pvalue_shapiro:.2f}")
if pvalue_shapiro > 0.05:
    print("Para Amostras controle: Não rejeitamos H0: Os dados são normais (p > 0.05)")
else:
    print("Para Amostras controle:Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)")
# Para Amostras controle: Não rejeitamos H0: Os dados são normais (p > 0.05)

stat_shapiro, pvalue_shapiro = shapiro(gdm['oGTT_1h'])
print(f"Stat_shapiro:{stat_shapiro:.2f}, p_value:{pvalue_shapiro:.2f}")
if pvalue_shapiro > 0.05:
    print("Para Amostras GDM :Não rejeitamos H0: Os dados são normais (p > 0.05)")
else:
    print("Para Amostras GDM :Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)")    
#Para Amostras GDM :Não rejeitamos H0: Os dados são normais (p > 0.05)


#Homocedastticidade
stat_levene, pvalue_levene = levene(control['oGTT_1h'], gdm['oGTT_1h'])
print(f"Stat_levene:{stat_shapiro:.2f}, p_value:{pvalue_levene:.2f}")
if pvalue_levene > 0.05:
    print("Não rejeitamos H0: As variâncias são homogêneas (p > 0.05)")
else:
    print("Rejeitamos H0: As variâncias NÃO são homogêneas (p ≤ 0.05)")
#Rejeitamos H0: As variâncias NÃO são homogêneas (p ≤ 0.05)    

#Teste t de Welch (correction =False)
welch_test_ogtt = pg.ttest(control['oGTT_1h'], gdm['oGTT_1h'],
                           correction=False)
print(welch_test_ogtt)

#Estatística do teste: -4.470, p-valor: 0.00010
#Há diferença significativa entre os grupos (p ≤ 0.05)

print('-----------------------------------------------------')

#### PESO
#NORMALIDADE
stat_shapiro, pvalue_shapiro = shapiro(control['Peso'])
print(f"Stat_shapiro:{stat_shapiro:.2f}, p_value:{pvalue_shapiro:.2f}")
if pvalue_shapiro > 0.05:
    print("Para Amostras controle: Não rejeitamos H0: Os dados são normais (p > 0.05)")
else:
    print("Para Amostras controle:Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)")
#Para Amostras controle: Não rejeitamos H0: Os dados são normais (p > 0.05)

stat_shapiro, pvalue_shapiro = shapiro(gdm['Peso'])
print(f"Stat_shapiro:{stat_shapiro:.2f}, p_value:{pvalue_shapiro:.2f}")
if pvalue_shapiro > 0.05:
    print("Para Amostras GDM :Não rejeitamos H0: Os dados são normais (p > 0.05)")
else:
    print("Para Amostras GDM :Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)")   
#Para Amostras GDM :Não rejeitamos H0: Os dados são normais (p > 0.05)

#Homocedastticidade
stat_levene, pvalue_levene = levene(control['Peso'],
                                    gdm['Peso'])
print(f"Stat_levene:{stat_shapiro:.2f}, p_value:{pvalue_levene:.2f}")
if pvalue_levene > 0.05:
    print("Não rejeitamos H0: As variâncias são homogêneas (p > 0.05)")
else:
    print("Rejeitamos H0: As variâncias NÃO são homogêneas (p ≤ 0.05)")
#Rejeitamos H0: As variâncias NÃO são homogêneas (p ≤ 0.05)    

# Teste t de Welch (correction =False)
welch_test_peso = pg.ttest(control['Peso'], gdm['Peso'],
                           correction=False) 
print(welch_test_peso)
#Estatística do teste: -3.114, p-valor: 0.00395
#Há diferença significativa entre os grupos (p ≤ 0.05)    

print('-----------------------------------------------------')
  
#### CSH genes 
#NORMALIDADE
stat_shapiro, pvalue_shapiro = shapiro(control['CSH_group'])
print(f"Stat_shapiro:{stat_shapiro:.2f}, p_value:{pvalue_shapiro:.2f}")
if pvalue_shapiro > 0.05:
    print("Para Amostras controle: Não rejeitamos H0: Os dados são normais (p > 0.05)")
else:
    print("Para Amostras controle:Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)")
#Para Amostras controle:Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05) 

stat_shapiro, pvalue_shapiro = shapiro(gdm['CSH_group'])
print(f"Stat_shapiro:{stat_shapiro:.2f}, p_value:{pvalue_shapiro:.2f}")
if pvalue_shapiro > 0.05:
    print("Para Amostras GDM :Não rejeitamos H0: Os dados são normais (p > 0.05)")
else:
    print("Para Amostras GDM :Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)")    
#Para Amostras GDM :Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)

#Homocedastticidade
stat_levene, pvalue_levene = levene(control['CSH_group'],
                                    gdm['CSH_group'])
print(f"Stat_levene:{stat_shapiro:.2f}, p_value:{pvalue_levene:.2f}")
if pvalue_levene > 0.05:
    print("Não rejeitamos H0: As variâncias são homogêneas (p > 0.05)")
else:
    print("Rejeitamos H0: As variâncias NÃO são homogêneas (p ≤ 0.05)")
#Não rejeitamos H0: As variâncias são homogêneas (p > 0.05)   

# Teste Mann Whiteney
mwu_csh = pg.mwu(control['CSH_group'], gdm['CSH_group'], alternative='two-sided')
print(mwu_csh)

#Estatística U: 143.0, p-valor: 0.0828
#Não há diferença significativa entre os grupos (p > 0.05)

print('-----------------------------------------------------')

#### IMC
#NORMALIDADE
stat_shapiro, pvalue_shapiro = shapiro(control['IMC'])
print(f"Stat_shapiro:{stat_shapiro:.2f}, p_value:{pvalue_shapiro:.2f}")
if pvalue_shapiro > 0.05:
    print("Para Amostras controle: Não rejeitamos H0: Os dados são normais (p > 0.05)")
else:
    print("Para Amostras controle:Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)")
#Para Amostras controle: Não rejeitamos H0: Os dados são normais (p > 0.05)

stat_shapiro, pvalue_shapiro = shapiro(gdm['IMC'])
print(f"Stat_shapiro:{stat_shapiro:.2f}, p_value:{pvalue_shapiro:.2f}")
if pvalue_shapiro > 0.05:
    print("Para Amostras GDM :Não rejeitamos H0: Os dados são normais (p > 0.05)")
else:
    print("Para Amostras GDM :Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)")   
#Para Amostras GDM :Não rejeitamos H0: Os dados são normais (p > 0.05)

#Homocedastticidade
stat_levene, pvalue_levene = levene(control['IMC'],
                                    gdm['IMC'])
print(f"Stat_levene:{stat_shapiro:.2f}, p_value:{pvalue_levene:.2f}")
if pvalue_levene > 0.05:
    print("Não rejeitamos H0: As variâncias são homogêneas (p > 0.05)")
else:
    print("Rejeitamos H0: As variâncias NÃO são homogêneas (p ≤ 0.05)")
#Não rejeitamos H0: As variâncias são homogêneas (p > 0.05)
      
# Teste t de Student (equal_var=True)

t_test_IMC = pg.ttest(control['IMC'], gdm['IMC'], correction=True) 
print(t_test_IMC)
#Estatística do teste: -3.138, p-valor: 0.00323
#Há diferença significativa entre os grupos (p ≤ 0.05)

print('-----------------------------------------------------')

#### oGTT_2h
#NORMALIDADE
stat_shapiro, pvalue_shapiro = shapiro(control2['ogtt_2h'])
print(f"Stat_shapiro:{stat_shapiro:.2f}, p_value:{pvalue_shapiro:.2f}")
if pvalue_shapiro > 0.05:
    print("Para Amostras controle: Não rejeitamos H0: Os dados são normais (p > 0.05)")
else:
    print("Para Amostras controle:Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)")
#Para Amostras controle: Não rejeitamos H0: Os dados são normais (p > 0.05)

stat_shapiro, pvalue_shapiro = shapiro(gdm2['ogtt_2h'])
print(f"Stat_shapiro:{stat_shapiro:.2f}, p_value:{pvalue_shapiro:.2f}")
if pvalue_shapiro > 0.05:
    print("Para Amostras GDM :Não rejeitamos H0: Os dados são normais (p > 0.05)")
else:
    print("Para Amostras GDM :Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)")    
#Para Amostras GDM :Não rejeitamos H0: Os dados são normais (p > 0.05)

#Homocedastticidade
stat_levene, pvalue_levene = levene(control2['ogtt_2h'],
                                    gdm2['ogtt_2h'])
print(f"Stat_levene:{stat_shapiro:.2f}, p_value:{pvalue_levene:.2f}")
if pvalue_levene > 0.05:
    print("Não rejeitamos H0: As variâncias são homogêneas (p > 0.05)")
else:
    print("Rejeitamos H0: As variâncias NÃO são homogêneas (p ≤ 0.05)")
#Não rejeitamos H0: As variâncias são homogêneas (p > 0.05)    

#Teste t de Student (equal_var=True)
t_test_ogtt2 = pg.ttest(control2['ogtt_2h'], gdm2['ogtt_2h'], correction=True)

#Estatística do teste: -4.016, p-valor: 0.00026
#Há diferença significativa entre os grupos (p ≤ 0.05)
    
print('-----------------------------------------------------')


#### Altura
#NORMALIDADE
stat_shapiro, pvalue_shapiro = shapiro(control2['maternal_height'])
print(f"Stat_shapiro:{stat_shapiro:.2f}, p_value:{pvalue_shapiro:.2f}")
if pvalue_shapiro > 0.05:
    print("Para Amostras controle: Não rejeitamos H0: Os dados são normais (p > 0.05)")
else:
    print("Para Amostras controle:Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)")
#Para Amostras controle: Não rejeitamos H0: Os dados são normais (p > 0.05)  

stat_shapiro, pvalue_shapiro = shapiro(gdm2['maternal_height'])
print(f"Stat_shapiro:{stat_shapiro:.2f}, p_value:{pvalue_shapiro:.2f}")
if pvalue_shapiro > 0.05:
    print("Para Amostras GDM :Não rejeitamos H0: Os dados são normais (p > 0.05)")
else:
    print("Para Amostras GDM :Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)") 
#Para Amostras GDM :Não rejeitamos H0: Os dados são normais (p > 0.05)

#Homocedastticidade
stat_levene, pvalue_levene = levene(control2['maternal_height'],
                                    gdm2['maternal_height'])
print(f"Stat_levene:{stat_shapiro:.2f}, p_value:{pvalue_levene:.2f}")
if pvalue_levene > 0.05:
    print("Não rejeitamos H0: As variâncias são homogêneas (p > 0.05)")
else:
    print("Rejeitamos H0: As variâncias NÃO são homogêneas (p ≤ 0.05)")
#Não rejeitamos H0: As variâncias são homogêneas (p > 0.05)
    
# Teste t de Student (equal_var=True)
t_test_altura = pg.ttest(control2['maternal_height'],
                         gdm2['maternal_height'], correction=True)

#Estatística do teste: -0.592, p-valor: 0.55708
#Não há diferença significativa entre os grupos (p > 0.05)

print('-----------------------------------------------------')
    
#### Idade 
#NORMALIDADE
stat_shapiro, pvalue_shapiro = shapiro(control2['maternal_age'])
print(f"Stat_shapiro:{stat_shapiro:.2f}, p_value:{pvalue_shapiro:.2f}")
if pvalue_shapiro > 0.05:
    print("Para Amostras controle: Não rejeitamos H0: Os dados são normais (p > 0.05)")
else:
    print("Para Amostras controle:Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)")
#Para Amostras controle: Não rejeitamos H0: Os dados são normais (p > 0.05)

stat_shapiro, pvalue_shapiro = shapiro(gdm2['maternal_age'])
print(f"Stat_shapiro:{stat_shapiro:.2f}, p_value:{pvalue_shapiro:.2f}")
if pvalue_shapiro > 0.05:
    print("Para Amostras GDM :Não rejeitamos H0: Os dados são normais (p > 0.05)")
else:
    print("Para Amostras GDM :Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)")   
#Para Amostras GDM :Rejeitamos H0: Os dados NÃO são normais (p ≤ 0.05)

#Homocedasticidade
stat_levene, pvalue_levene = levene(control2['maternal_age'],
                                    gdm2['maternal_age'])
print(f"Stat_levene:{stat_shapiro:.2f}, p_value:{pvalue_levene:.2f}")
if pvalue_levene > 0.05:
    print("Não rejeitamos H0: As variâncias são homogêneas (p > 0.05)")
else:
    print("Rejeitamos H0: As variâncias NÃO são homogêneas (p ≤ 0.05)")
#Não rejeitamos H0: As variâncias são homogêneas (p > 0.05)
    
# Teste Mann-Whitney

mwu_csh_idade = pg.mwu(control2['maternal_age'],
                                    gdm2['maternal_age'], alternative='two-sided')


#Estatística U: 192.0, p-valor: 0.6460
#Não há diferença significativa entre os grupos (p > 0.05)  



print('-----------------------------------------------------')

#### IMC Categórico
# Tabela de Contingência
contingency_table = pd.crosstab(data_stats['group'], data_stats['BMI_adjust_cat'])
print(contingency_table)

# Teste exato de Fisher 
odds_ratio, p_value_fisher = fisher_exact(contingency_table)
print("OR:", odds_ratio, "p:", p_value_fisher)

res = sm.stats.Table2x2(contingency_table.values)
ci_low, ci_upp = res.oddsratio_confint(alpha=0.05)
print("IC 95% do OR:", (ci_low, ci_upp))

#Diferença estatisticamente significativa (p = 0.0036)
#Pacientes com GDM têm 17.3 vezes mais chances de estarem na categoria Sobrepeso

#%% Gráficos 

palette_order = {"Control": "lightblue", "GDM": "salmon"}

#### oGTT_1h
plt.figure(figsize =(4,5), dpi=300)
sns.boxplot(data_stats,
            x='group',
            y = 'oGTT_1h',
            order = ['Control', 'GDM'],
            palette=palette_order,
            showmeans = True,
            meanprops={"marker": "x",
                       "markerfacecolor": "black",
                       "markeredgecolor": "black"},
            boxprops=dict(edgecolor="black"),  
            whiskerprops=dict(color="black"),  
            capprops=dict(color="black"),  
            medianprops=dict(color="black")) 
plt.title('1h post-OGTT',
             fontsize=18,
             fontname='Arial',
             fontweight='bold',
             y=1.05,
             pad =20)
plt.ylabel('Glycemia\n(mmol/L)',
           fontsize = 16,
           fontweight ='bold',
           fontname = 'Arial')
plt.xlabel('')
plt.xticks([0,1] , labels = ['Control','GDM'],
          fontsize = 16,
          fontname = 'Arial',
          fontweight ='bold')
plt.ylim(0 , 16)
plt.yticks(fontsize = 16,
           fontname = 'Arial')

# Adicionar linha de significância
y_max = data_stats['oGTT_1h'].max()  
plt.plot([0, 1], [y_max + 1, y_max + 1],
         color='black', lw=1)  
plt.text(0.5, y_max + 1.25, 'p=0.0001',
         ha='center', va='bottom',
         color='black', fontsize=12)  

sns.despine()
plt.tight_layout()
plt.show()

#### Peso
plt.figure(figsize =(4,5), dpi=300)
sns.boxplot(data_stats,
            x='group',
            y = 'Peso',
            order = ['Control', 'GDM'],
            palette=palette_order,
            showmeans = True,
            meanprops={"marker": "x",
                       "markerfacecolor": "black",
                       "markeredgecolor": "black"},
            boxprops=dict(edgecolor="black"),  
            whiskerprops=dict(color="black"),  
            capprops=dict(color="black"),  
            medianprops=dict(color="black")) 
plt.title('Peso',
             fontsize=18,
             fontname='Arial',
             fontweight='bold',
             y=1.05,
             pad =20)
plt.ylabel('Kilogramas (Kg)',
           fontsize = 16,
           fontweight ='bold',
           fontname = 'Arial')
plt.xlabel('')
plt.xticks([0,1] , labels = ['Controle','DMG'],
          fontsize = 16,
          fontname = 'Arial',
          fontweight ='bold')
plt.ylim(0 , 100)
plt.yticks(fontsize = 16,
           fontname = 'Arial')

# Adicionar linha de significância
y_max = data_stats['Peso'].max()  
plt.plot([0, 1], [y_max + 5, y_max + 5],
         color='black', lw=1)  
plt.text(0.5, y_max + 5, 'p= 0.004',
         ha='center', va='bottom',
         color='black', fontsize=12)  

sns.despine()
plt.tight_layout()
plt.show()

#### CSH
plt.figure(figsize =(4,5), dpi=300)
sns.boxplot(data_stats,
            x='group',
            y = 'CSH_group',
            order = ['Control', 'GDM'],
            palette=palette_order,
            showmeans = True,
            meanprops={"marker": "x",
                       "markerfacecolor": "black",
                       "markeredgecolor": "black"},
            boxprops=dict(edgecolor="black"),  
            whiskerprops=dict(color="black"),  
            capprops=dict(color="black"),  
            medianprops=dict(color="black")) 
plt.title('Grupo CSH de Genes',
             fontsize=18,
             fontname='Arial',
             fontweight='bold',
             y=1.05,
             pad =20)
plt.ylabel('Z-Score de Log(CPM + 1)',
           fontsize = 16,
           fontweight ='bold',
           fontname = 'Arial')
plt.xlabel('')
plt.xticks([0,1] , labels = ['Controle','DMG'],
          fontsize = 16,
          fontname = 'Arial',
          fontweight ='bold')
plt.ylim(-8 , 3)
plt.yticks(fontsize = 16,
           fontname = 'Arial')

# Adicionar linha de significância
y_max = data_stats['CSH_group'].max()  
plt.plot([0, 1], [y_max + 0.20, y_max + 0.20],
         color='black', lw=1)  
plt.text(0.5, y_max + 0.22, 'p=0.08',
         ha='center', va='bottom',
         color='black', fontsize=12)  
sns.despine()
plt.tight_layout()
plt.show()

#### OGTT2 
plt.figure(figsize =(4,5), dpi=300)
sns.boxplot(data_2,
            x='group',
            y = 'ogtt_2h',
            order = ['Control', 'GDM'],
            palette=palette_order,
            showmeans = True,
            meanprops={"marker": "x",
                       "markerfacecolor": "black",
                       "markeredgecolor": "black"},
            boxprops=dict(edgecolor="black"),  
            whiskerprops=dict(color="black"),  
            capprops=dict(color="black"),  
            medianprops=dict(color="black")) 
plt.title('TOTG-2h',
             fontsize=18,
             fontname='Arial',
             fontweight='bold',
             y=1.05,
             pad =20)
plt.ylabel('Glicemia sanguínea\n(mmol/L)',
           fontsize = 16,
           fontweight ='bold',
           fontname = 'Arial')
plt.xlabel('')
plt.xticks([0,1] , labels = ['Controle','DMG'],
          fontsize = 16,
          fontname = 'Arial',
          fontweight ='bold')
plt.ylim(0 , 16)
plt.yticks(fontsize = 16,
           fontname = 'Arial')
# Adicionar linha de significância
y_max = data_2['ogtt_2h'].max()  
plt.plot([0, 1], [y_max + 1, y_max + 1],
         color='black', lw=1)  
plt.text(0.5, y_max + 1.2, 'p=0.0002',
         ha='center', va='bottom',
         color='black', fontsize=12)  
sns.despine()
plt.tight_layout()
plt.show()

#### Altura
plt.figure(figsize =(4,5), dpi=300)
sns.boxplot(data_2,
            x='group',
            y = 'maternal_height',
            order = ['Control', 'GDM'],
            palette=palette_order,
            showmeans = True,
            meanprops={"marker": "x",
                       "markerfacecolor": "black",
                       "markeredgecolor": "black"},
            boxprops=dict(edgecolor="black"),  
            whiskerprops=dict(color="black"),  
            capprops=dict(color="black"),  
            medianprops=dict(color="black")) 
plt.title('Altura',
             fontsize=18,
             fontname='Arial',
             fontweight='bold',
             y=1.05,
             pad =20)
plt.ylabel('Centímetros (cm)',
           fontsize = 16,
           fontweight ='bold',
           fontname = 'Arial')
plt.xlabel('')
plt.xticks([0,1] , labels = ['Controle','DMG'],
          fontsize = 16,
          fontname = 'Arial',
          fontweight ='bold')
plt.ylim(140 , 180)
plt.yticks(fontsize = 16,
           fontname = 'Arial')
# Adicionar linha de significância
y_max = data_2['maternal_height'].max()  
plt.plot([0, 1], [y_max + 2, y_max + 2],
         color='black', lw=1)  
plt.text(0.5, y_max + 2.5, 'p=0.557',
         ha='center', va='bottom',
         color='black', fontsize=12)  
sns.despine()
plt.tight_layout()
plt.show() 

#### Idade
plt.figure(figsize =(4,5), dpi=300)
sns.boxplot(data_2,
            x='group',
            y = 'maternal_age',
            order = ['Control', 'GDM'],
            palette=palette_order,
            showmeans = True,
            meanprops={"marker": "x",
                       "markerfacecolor": "black",
                       "markeredgecolor": "black"},
            boxprops=dict(edgecolor="black"),  
            whiskerprops=dict(color="black"),  
            capprops=dict(color="black"),  
            medianprops=dict(color="black")) 
plt.title('Idade',
             fontsize=18,
             fontname='Arial',
             fontweight='bold',
             y=1.05,
             pad =20)
plt.ylabel('Anos',
           fontsize = 16,
           fontweight ='bold',
           fontname = 'Arial')
plt.xlabel('')
plt.xticks([0,1] , labels = ['Controle','DMG'],
          fontsize = 16,
          fontname = 'Arial',
          fontweight ='bold')
plt.ylim(22 , 46)
plt.yticks(np.arange(22, 48, 4), 
           fontsize=16, fontname='Arial')
# Adicionar linha de significância
y_max = data_2['maternal_age'].max()  
plt.plot([0, 1], [y_max + 1, y_max + 1],
         color='black', lw=1)  
plt.text(0.5, y_max + 1.2, 'p=0.646',
         ha='center', va='bottom',
         color='black', fontsize=12)  
sns.despine()
plt.tight_layout()
plt.show()

#### IMC
plt.figure(figsize =(4,5), dpi=300)
sns.boxplot(data,
            x='group',
            y = 'IMC',
            order = ['Control', 'GDM'],
            palette=palette_order,
            showmeans = True,
            meanprops={"marker": "x",
                       "markerfacecolor": "black",
                       "markeredgecolor": "black"},
            boxprops=dict(edgecolor="black"),  
            whiskerprops=dict(color="black"),  
            capprops=dict(color="black"),  
            medianprops=dict(color="black")) 
plt.title('Body Mass Index',
             fontsize=18,
             fontname='Arial',
             fontweight='bold',
             y=1.05,
             pad =20)
plt.ylabel('kg/m2',
           fontsize = 16,
           fontweight ='bold',
           fontname = 'Arial')
plt.xlabel('')
plt.xticks([0,1] , labels = ['Control','GDM'],
          fontsize = 16,
          fontname = 'Arial',
          fontweight ='bold')
plt.ylim(0 , 35)
plt.yticks(fontsize = 16,
           fontname = 'Arial')
# Adicionar linha de significância
y_max = data['IMC'].max()  
plt.plot([0, 1], [y_max + 1, y_max + 1],
         color='black', lw=1)  
plt.text(0.5, y_max + 1.2, 'p=0.003',
         ha='center', va='bottom',
         color='black', fontsize=12)  
sns.despine()
plt.tight_layout()
plt.show()

#### Categorias
palette = {"Control": '#a6cee3', "GDM": '#fdbf6f'}
colors = [palette['Control'], palette['GDM']]  
edge_color = 'black'  

plt.figure(figsize=(4,4), dpi=300)
cross_tab = pd.crosstab(data['group'], 
                       data['BMI_adjust_cat'], 
                       normalize='index')
ax = cross_tab.plot(kind='bar', 
                   stacked=True, 
                   color=['#a6cee3','#fdbf6f'],  
                   edgecolor=edge_color,
                   linewidth=1)  
plt.title('BMI Categories',
         fontsize=16,
         fontname='Arial',
         fontweight='bold',
         y=1.05,
         pad=20)
plt.ylabel('Proportion',
          fontsize=14,
          fontweight='bold',
          fontname='Arial')
plt.xlabel('',
          fontsize=14,
          fontweight='bold',
          fontname='Arial')
plt.xticks([0, 1], 
           labels=['Control', 'GDM'],
           fontsize=14,
           fontname='Arial',
           fontweight='bold',
           rotation=0)
plt.yticks(fontsize=14,
          fontname='Arial')
plt.legend(title='',
           labels =['Normal weight', 'Overweight'],
          bbox_to_anchor=(1,1),
          frameon=False,
          fontsize=12,
          title_fontsize=14)
sns.despine()
plt.tight_layout()
plt.show()
