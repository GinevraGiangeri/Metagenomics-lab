#import libraries
import numpy as np
import pandas as pd
import os, glob, fnmatch, collections
import matplotlib.pyplot as plt

plt.style.use("ggplot")

#let the user choose between KEGG ortologues or modules and the taxonomic level to use in the stratification step
anno_selection = str(input("Choose between KEGG orthologues (ko) or KEGG modules (mo) : "))
print("\n")
taxa_selection = str(input("Choose the taxonomic level among Domain, Phylum, Class, Order, Family or Species : "))
print("\n")

#import MAGs abundances
path1 = "/home/gabriele/Documents/Lab_metagenomics/Data1/MAGs_coverage.txt"
ab = pd.read_csv(path1, sep='\t', engine='python')

df1 = ab['Bin Id']
df2 = ab.iloc[:,ab.columns.str.contains('populations')]
abundances = pd.concat([df1,df2], axis = 1)
mean_ab = abundances.sort_values(by='Bin Id')

mean_ab.to_csv(r'/home/gabriele/Documents/Lab_metagenomics/Test/MAGs_abundances.csv', sep='\t')

#import taxonomic assignment
path2 = "/home/gabriele/Documents/Lab_metagenomics/Data1/taxonomy.txt"
taxonomy = pd.read_csv(path2, sep='\t', engine='python')

new = taxonomy["NCBI classification"].str.split(";", n = 7, expand = True) 
taxonomy["Domain"]= new[0]
taxonomy["Phylum"]= new[1]
taxonomy["Class"]= new[2]
taxonomy["Order"]= new[3]
taxonomy["Family"]= new[4]
taxonomy["Genus"]= new[5]
taxonomy["Species"]= new[6]
taxonomy.drop(columns =["NCBI classification"], inplace = True)

taxonomy['Class'].loc[(taxonomy['Class'] == 'c__')] = 'Unclassified'
taxonomy['Order'].loc[(taxonomy['Order'] == 'o__')] = 'Unclassified'
taxonomy['Family'].loc[(taxonomy['Family'] == 'f__')] = 'Unclassified'
taxonomy['Genus'].loc[(taxonomy['Genus'] == 'g__')] = 'Unclassified'
taxonomy['Species'].loc[(taxonomy['Species'] == 's__')] = 'Unclassified'

taxonomy = taxonomy.sort_values(by='user_genome')
taxonomy = taxonomy.filter(['user_genome','Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'], axis=1)
taxonomy.to_csv(r'/home/gabriele/Documents/Lab_metagenomics/Test/Tassonomic_assignment.csv', sep='\t')

# import eggNOG annotation files and calculate the total weighted occurrences for all the MAGs annotation files
path = "/home/gabriele/Documents/Lab_metagenomics/Data1/MAGs_annotation"
os.chdir(path)
cwd = os.getcwd()
csv_files = glob.glob(os.path.join(cwd, "*.annotations")) #group together the MAG annotation files
list_ann = list()
column_list = mean_ab.columns.tolist()
column_list.remove('Bin Id')

for element in column_list:
    KEGG_list = list()
    weight_list = list()
    weight_list = mean_ab[element].tolist()
    species_list = list()
    species_list = mean_ab['Bin Id'].tolist()
    db_list = pd.DataFrame()

    for f in csv_files:
            list_ann.append(f)
            df = pd.read_csv(f, sep='\t', header=3)
            df.drop(df.tail(3).index, inplace = True)
            df = df.loc[:, ['KEGG_ko', 'KEGG_Module']] 
            df = df.astype(str)
            
            if anno_selection == 'ko':
                df = df['KEGG_ko'].str.split(',', expand=True).stack().reset_index(drop=True, level=1)
            elif anno_selection == 'mo':
                df = df['KEGG_Module'].str.split(',', expand=True).stack().reset_index(drop=True, level=1)
            else:
                print("\nERROR: Please type 'ko' if you want to select the kegg ortologues or 'mo' if you want the kegg module")
                break
            
            KEGG_list = df.tolist() #list of all the kegg codes
            counts = list(collections.Counter(KEGG_list).items()) #create a list of tuples: kegg_code , occurrences
            
            d1 = dict()
            for KEGG,n in counts:
                d1.setdefault(KEGG, []).append(n) #create a dictionary from the list of tuples 'counts' 
    
            kegg = pd.DataFrame.from_dict(d1, orient='index', columns=['Occurrences']) #create the dataframe from the dictionary
            kegg = kegg.drop('nan') #remove the kegg module not detected, identified as NaN
            index_names = kegg[ kegg['Occurrences'] == 1 ].index #remove all the kegg module with only 1 occurrence
            kegg.drop(index_names, inplace = True) 
            kegg.reset_index(level=0, inplace=True)
            
            index = list_ann.index(f)
            kegg['Weight'] = weight_list[index] #add a column with the weight for that species    
            kegg['Occurrences_weighted'] = kegg['Occurrences'].astype(int)*kegg['Weight'].astype(float) #weight the occurrences
            
            kegg = kegg.sort_values(by = 'Occurrences_weighted', ascending=False) #sort in descending order
            kegg = kegg.filter(['index', 'Occurrences_weighted']) #keep only the columns we need
            
            decimals = 1    
            kegg['Occurrences_weighted'] = kegg['Occurrences_weighted'].apply(lambda x: round(x, decimals)) #round

            db_list = db_list.append(kegg, sort=False) #create a unique databases

    db = db_list.groupby('index')['Occurrences_weighted'].sum().sort_values(ascending=False) #sum and sort occurrences

    plt.figure(figsize=(10, 8))
    db.head(10).plot.bar(color="darkseagreen")
    plt.title("KEGG codes ordered by total occurrences {} (Top10)".format(element))
    plt.xlabel("Number of occurrences")
    plt.ylabel("KEGG codes")
    plt.savefig('/home/gabriele/Documents/Lab_metagenomics/Test/{}.png'.format(element), dpi=300, bbox_inches='tight')

    print("Weighted abundances of {} finished".format(element))

#stratification
d = {}
d_taxa = {}
print("\nStratification is running ...")
for f in csv_files:
    index = list_ann.index(f)
    s = species_list[index]
    df = pd.read_csv(f, sep='\t', header=3)
    df.drop(df.tail(3).index, inplace = True)
    df = df.loc[:, ['KEGG_ko', 'KEGG_Module']] 
    df = df.astype(str)
        
    if anno_selection == 'ko':
        df = df['KEGG_ko'].str.split(',', expand=True).stack().reset_index(drop=True, level=1)
    elif anno_selection == 'mo':
        df = df['KEGG_Module'].str.split(',', expand=True).stack().reset_index(drop=True, level=1)
    else:
        print("\nERROR: Please type 'ko' if you want to select the kegg ortologues or 'mo' if you want the kegg module")
        break
        
    KEGG_list = df.tolist()
    
    for code in KEGG_list:
        if (code not in d.keys() and code!="nan"):
            l= taxonomy[taxonomy['user_genome'] == s.replace('.emapper.annotations', '')][taxa_selection].values
            l = str(l).replace("['","").replace("']","")
            d[code] = [l]
        elif (code in d.keys() and code!="nan"):
            l = taxonomy[taxonomy['user_genome'] == s.replace('.emapper.annotations', '')][taxa_selection].values
            l = str(l).replace("['","").replace("']","")
            n = d[code].append(l)
        else:
            continue

taxa_list = taxonomy[taxa_selection].tolist()
x = collections.Counter(taxa_list).items()

for taxa,n in x:
    d_taxa.setdefault(taxa, []).append(n)

taxa = pd.DataFrame.from_dict(d_taxa, orient='index', columns=['Occurrences'])
taxa['Percentages'] = (taxa['Occurrences']/taxa['Occurrences'].sum())*100
taxa.reset_index(level=0, inplace=True)
taxa = taxa.sort_values(by='index')

for i in d.keys():
    list_names = d[i]
    frequency = collections.Counter(list_names).items()
    data = pd.DataFrame(frequency, columns =[taxa_selection, 'Frequency'])
    data = data.sort_values(by=taxa_selection)
    for j in range(len(data)):
        data['Real Frequency'] = data['Frequency'].astype(float)*taxa['Percentages'].astype(float)

    data.to_csv(r'/home/gabriele/Documents/Lab_metagenomics/Test/Stratification/{}.csv'.format(i), sep='\t', index=False, float_format='%.2f')
    
print("\nSTRIM run succefully ! Check the output directory to inspect the results")