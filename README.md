# ORAL-MEDICATIONS-SARS
This is a small drug discovery bioinformatics project. We are searching for molecules that could be taken orally by a human to treat SARS.
It is my first complete protein delocalisation from chembl 
get the code on Github.
pip install chembl_webresource_client
%matplotlib inline
import matplotlib.pyplot as plt
import sys
import os
sys.path.append('/usr/local/lib/python3.7/site-packages/')  ### This is needed to make sure we can import RDKit in Google Colab.
import pandas as pd
from chembl_webresource_client.new_client import new_client
target = new_client.target
target_query = target.search('coronavirus')
targets = pd.DataFrame.from_dict(target_query)
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', None)
targets
selected_target = targets.target_chembl_id[4]
selected_target
activity = new_client.activity
res = activity.filter(target_chembl_id=selected_target)
df = pd.DataFrame.from_dict(res)
Narrowing down to IC50
activity2 = new_client.activity
res2 = activity2.filter(target_chembl_id=selected_target).filter(type='IC50')
df2 = pd.DataFrame.from_dict(res2)
df2.head(10)
Let's save a CSV at this point and the reload the new CSV.
df2.to_csv('IMP13_bioactivity_data.csv', index=False)
df3 = df2
### df3 = pd.read_csv(r'C:\Users\PheXeRiaN\Desktop\Drug Discovery\IMP13_bioactivity_data.csv')   ### This is command for when saving csv not in Google Colab to load from CSV.
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', None)
df3.head(5)
df4 = df3[df3.value.notna()]
df4.head(5)
bioactivity_class = []
for x in df4.standard_value:
    if float(x) >= 10000:
        bioactivity_class.append('inactive')
    elif float(x) <= 1000:
        bioactivity_class.append('active')
    else:
        bioactivity_class.append('intermediate')
Now we need to remove duplicates for a few of the rows and moves features around to a new dataframe.
[ ]
molecule_chembl_id = []
for x in df4.molecule_chembl_id:
    molecule_chembl_id.append(x)
[ ]
canonical_smiles = []
for x in df4.canonical_smiles:
    canonical_smiles.append(x)

[ ]
standard_value = []
for x in df4.standard_value:
    standard_value.append(x)
[ ]
data_tuples = list(zip(molecule_chembl_id, canonical_smiles, bioactivity_class, standard_value))

df5 = pd.DataFrame(data_tuples, columns=['molecule_chembl_id', 'canonical_smiles', 'bioactivity_class', 'standard_value'])
df5.head(5)
[ ]
df5.to_csv('bioactivity_preprocessed_data.csv', index=False)

Calculate Lipinski descriptors. Rule of Five for ADME pharmacokinetic profile.

! wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.8.2-Linux-x86_64.sh
! chmod +x Miniconda3-py37_4.8.2-Linux-x86_64.sh
! bash ./Miniconda3-py37_4.8.2-Linux-x86_64.sh -b -f -p /usr/local
! conda install -c rdkit rdkit -y
import sys
sys.path.append('/usr/local/lib/python3.7/site-packages/')
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski

def lipinski(smiles, verbose=False):

    moldata= []
    for elem in smiles:
        mol=Chem.MolFromSmiles(elem) 
        moldata.append(mol)
       
    baseData= np.arange(1,1)
    i=0  
    for mol in moldata:        
       
        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)
           
        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])   
    
        if(i==0):
            baseData=row
        else:
            baseData=np.vstack([baseData, row])
        i=i+1      
    
    columnNames=["MW","LogP","NumHDonors","NumHAcceptors"]   
    descriptors = pd.DataFrame(data=baseData,columns=columnNames)
    
    return descriptors

### Made from code from https://codeocean.com/capsule/8848590/tree/v1

df_lipinski = lipinski(df.canonical_smiles)





[ ]


df_lipinski.head(10)

Combining dataframes

[ ]


df6 = pd.concat([df5, df_lipinski], axis=1)

[ ]


df6.head(10)

import seaborn as sns

plt.figure(figsize=(8,5))
sns.distplot(df6.standard_value)

Now need to convert IC50 to negative log scale (-log10(IC50)) which gives us a more uniform distribution for prediction.



[ ]


def pIC50(input):
    pIC50 = []

    for i in input['standard_value_norm']:
        molar = i*(10**-9) # Converts nanoMolar to Molar
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)
        
    return x







We will max out values at 100000000 to make it easier for us. This gets rid of negative values if negative log is calculated.

[ ]


df6.standard_value.dropna(inplace=True)





[ ]


df6.standard_value = pd.to_numeric(df6.standard_value).astype(float)





[ ]


def norm_value(input):
    norm = []

    for i in input['standard_value']:
        if i > 100000000:
          i = 100000000
        norm.append(i)

    input['standard_value_norm'] = norm
    x = input.drop('standard_value', 1)
        
    return x





[ ]


df_normalized = norm_value(df6)
df_normalized.head(10)







[ ]


df7 = pIC50(df_normalized)
df7.head(10)







[ ]


plt.figure(figsize=(8,5))
sns.distplot(df7.pIC50)







[ ]


df7.pIC50.describe()



count    105.000000 mean       4.876188 std        0.919190 min        3.000000 25%        4.346787 50%        4.823909 75%        5.055517 max        7.301030 Name: pIC50, dtype: float64






Removing intermediate class



[ ]


df_2class = df7[df7.bioactivity_class != 'intermediate']







Looking at Active vs Inactive molecules



[ ]


sns.set(style='darkgrid')
sns.countplot(x='bioactivity_class', data = df_2class, edgecolor='black')

plt.xlabel('Bioactivity class', fontsize = 17)
plt.ylabel('Frequency', fontsize = 17)
plt.title('Inactive vs Active Bar Graph', size = 20, fontweight = 'bold')

plt.savefig('plot_bioactivity_class.png')







[ ]


sns.scatterplot(x='MW', y='LogP', data=df_2class, hue='bioactivity_class', size='pIC50', edgecolor='black', alpha=0.7)

plt.xlabel('MW', fontsize=14)
plt.ylabel('LogP', fontsize=14)
plt.title('MW to LogP', size = 20, fontweight = 'bold')
plt.legend(bbox_to_anchor=(1, 1), loc=2, borderaxespad=0, fancybox = True)

plt.savefig('plot_MW_vs_LogP.png')







[ ]


sns.boxplot(x = 'bioactivity_class', y = 'pIC50', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14)
plt.ylabel('pIC50 value', fontsize=14)
plt.title('Boxplot of Inactive to Active', size = 20, fontweight = 'bold')

plt.savefig('plot_ic50.png')









Now a Mann-Whitney U Test to make sure both groups are statistically different. This MannWhitney test is for active vs inactive bioactivity.



[ ]


def mannwhitney(descriptor, verbose=True):
    from numpy.random import seed
    from numpy.random import randn
    from scipy.stats import mannwhitneyu


    seed(42)

### actives and inactives
    selection = [descriptor, 'bioactivity_class']
    df = df_2class[selection]
    active = df[df.bioactivity_class == 'active']
    active = active[descriptor]

    selection = [descriptor, 'bioactivity_class']
    df = df_2class[selection]
    inactive = df[df.bioactivity_class == 'inactive']
    inactive = inactive[descriptor]

### comparing samples
    stat, p = mannwhitneyu(active, inactive)
    print('Statistics=%.3f, p=%.3f' % (stat, p))

### interpret
    alpha = 0.05
    if p > alpha:
        interpretation = 'Same distribution (fail to reject H0)'
    else:
        interpretation = 'Different distribution (reject H0)'
  
    results = pd.DataFrame({'Descriptor':descriptor,
                          'Statistics':stat,
                          'p':p,
                          'alpha':alpha,
                          'Interpretation':interpretation}, index=[0])
    #filename = 'mannwhitneyu_' + descriptor + '.csv'
    #results.to_csv(filename)

    return results





[ ]


mannwhitney('pIC50')









Molecular Weight



[ ]


sns.boxplot(x = 'bioactivity_class', y = 'MW', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14)
plt.ylabel('MW', fontsize=14)
plt.title('MW to Bioactivity Class', size = 20, fontweight = 'bold')

plt.savefig('plot_MW.png')







[ ]


mannwhitney('MW')









LogP



[ ]


sns.boxplot(x = 'bioactivity_class', y = 'LogP', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14)
plt.ylabel('LogP', fontsize=14)
plt.title('LogP vs Bioactictivity Class', size = 20, fontweight='bold')

plt.savefig('plot_LogP.png')







[ ]


mannwhitney('LogP')









Number of Hydrogen Donors



[ ]


sns.boxplot(x = 'bioactivity_class', y = 'NumHDonors', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14)
plt.ylabel('LogP', fontsize=14)
plt.title('Hydrogen Donors vs Bioactictivity Class', size = 20, fontweight='bold')

plt.savefig('plot_NumHDonors.png')







[ ]


mannwhitney('NumHDonors')









Number of Hydrogen Acceptors



[ ]


sns.boxplot(x = 'bioactivity_class', y = 'NumHAcceptors', data = df_2class)

plt.xlabel('Bioactivity class', fontsize=14)
plt.ylabel('LogP', fontsize=14)
plt.title('Hydrogen Acceptors vs Bioactictivity Class', size = 20, fontweight='bold')

plt.savefig('plot_NumHacceptors.png')







[ ]


mannwhitney('NumHAcceptors')









Interpretations
pIC50 ------------ Statistically significant (This is expected as we pre-processed this data to split active and inactive)
MW -------------- Not significant
LogP ------------ Not significant
H Donors ------ Not significant
H Acceptors -- Not significant



[ ]


df_2class.to_csv('df_2class.csv', index=False)







We will now use PaDEL to calculate the molecular descriptors.



[ ]


! wget https://github.com/dataprofessor/bioinformatics/raw/master/padel.zip



df_2class.dropna(inplace=True)



/usr/local/lib/python3.6/dist-packages/ipykernel_launcher.py:1: SettingWithCopyWarning:  A value is trying to be set on a copy of a slice from a DataFrame  See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy   """Entry point for launching an IPython kernel. 




[ ]


selection = ['canonical_smiles','molecule_chembl_id']
df3_selection = df3[selection]
df3_selection.to_csv('molecule.smi', sep='\t', index=False, header=False)





[ ]


! cat molecule.smi | head -5



Cc1noc(C)c1CN1C(=O)C(=O)c2cc(C#N)ccc21	CHEMBL187579 O=C1C(=O)N(Cc2ccc(F)cc2Cl)c2ccc(I)cc21	CHEMBL188487 O=C1C(=O)N(CC2COc3ccccc3O2)c2ccc(I)cc21	CHEMBL185698 O=C1C(=O)N(Cc2cc3ccccc3s2)c2ccccc21	CHEMBL426082 O=C1C(=O)N(Cc2cc3ccccc3s2)c2c1cccc2[N+](=O)[O-]	CHEMBL187717 




[ ]


 cat molecule.smi | wc -l   ### We have 215 molecules to look transcribe.



105 






Calculating descriptors



[ ]


! cat padel.sh



java -Xms1G -Xmx1G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv 






Showing what is inside padel file in above line of code. This uses Java with 1 gig of memory. We removed the sodium and chloride (salts) which are in the chemical structure.



[ ]


! bash padel.sh




Preparing the X and Y Data Matrices
X Matrix



[ ]


df_X_matrix = pd.read_csv('descriptors_output.csv')





[ ]


df_X_matrix







[ ]


df_X_matrix = df_X_matrix.drop(columns=['Name'])
df_X_matrix









Y variable
COnverting IC50 to pIC50



[ ]


df_y = df7['pIC50']
df_y.dropna(inplace=True)





[ ]


df_y



0      5.142668 1      5.026872 2      4.869666 3      4.882397 4      5.698970 5      6.008774 6      5.316953 7      6.022276 8      4.950782 9      4.628932 10     4.900665 11     4.756962 12     4.346787 13     4.154902 14     4.180456 15     6.431798 16     4.903090 17     4.721246 18     4.602060 19     4.148742 20     5.958607 21     4.301030 22     5.522879 23     3.522879 24     3.602060 25     3.698970 26     4.000000 27     4.221849 28     4.346787 29     4.397940 30     4.823909 31     4.823909 32     4.920819 33     3.000000 34     3.301030 35     3.397940 36     3.455932 37     3.522879 38     3.522879 39     3.698970 40     3.698970 41     3.698970 42     3.698970 43     4.221849 44     4.397940 45     4.522879 46     4.823909 47     4.853872 48     4.958607 49     5.000000 50     6.045757 51     5.221849 52     4.920819 53     4.886057 54     4.886057 55     4.823909 56     4.795880 57     4.795880 58     4.795880 59     4.602060 60     4.494850 61     5.522879 62     5.301030 63     5.000000 64     4.823909 65     4.795880 66     4.744727 67     4.744727 68     4.698970 69     4.397940 70     6.522879 71     7.221849 72     7.200659 73     6.468521 74     6.568636 75     7.022276 76     7.187087 77     7.301030 78     6.769551 79     4.050122 80     4.605548 81     4.675718 82     3.644548 83     4.412289 84     4.841638 85     4.675718 86     5.795880 87     4.970616 88     5.036212 89     6.096910 90     5.055517 91     5.309804 92     4.522879 93     4.283997 94     4.057496 95     6.154902 96     5.920819 97     4.220404 98     4.767004 99     4.744727 100    4.974694 101    4.995679 102    4.939302 103    4.970616 104    4.102923 Name: pIC50, dtype: float64






Combining



[ ]


final_dataset = pd.concat([df_X_matrix,df_y], axis=1)





[ ]


missdata = df.isnull().values.any()
print(missdata)



True 




[ ]


if missdata == True:
    print('# of missing values:', final_dataset.isnull().values.sum())
else:
    print('No missing data')



# of missing values: 0 




[ ]


final_dataset.isnull().sum()



[ ]


final_dataset = final_dataset.dropna(axis=0)   ### Removing one final NaN value found in matrix.





[ ]


final_dataset







[ ]


final_dataset.to_csv('final_dataset,csv', index=False)







Modeling



[ ]


final_dataset = pd.read_csv('final_dataset,csv')





[ ]


final_dataset







[ ]


from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor





[ ]


X = final_dataset.drop('pIC50', axis=1).copy()
y = final_dataset.pIC50.copy()





[ ]


X.shape



(105, 881)




[ ]


y.shape



(105,)




[ ]


X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)





[ ]


model = RandomForestRegressor(n_estimators=100)
model.fit(X_train, y_train)
R_square = model.score(X_test, y_test)
R_square



0.532986696734584






R-squared values should be between 0.5 and 0.6 for most workable models. So this model is valid and usable.



[ ]


y_pred = model.predict(X_test)







Scatterplot of Predicted vs Experimental pIC50



[ ]


ax = sns.regplot(y_test, y_pred, scatter_kws = {'alpha':0.6})
ax.set_xlabel('Experimental pIC50', fontsize = 'large')
ax.set_ylabel('Predicted pIC50', fontsize='large')
plt.title('Predicted vs Experimental pIC50', size = 18, fontweight = 'bold')
ax.figure.set_size_inches(5, 5)
plt.show







[ ]


import plotly.express as px

fig = px.scatter(y_test, y_pred, color = "pIC50")

fig.update_layout(
    height=600,
    width=600,
    title_text='Predicted vs Experimental pIC50')

fig.show()







A good molecule to investigate would have a pIC50 of near 5.5 or greater. This would give a stronger inhibition of the binded protein on SARs virus.
This is the end. As you can see we can predict molecules with a pIC50 greater than 5.5, although a better goal to strive for is 6. Our model is valid with a r-squared value greater than 0.5 as well.
