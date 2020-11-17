import pandas as pd
import numpy as np
import re as re
# Standard library package
import io
import os
import requests as requests

# Import Biopython modules to interact with KEGG
from Bio import SeqIO
from Bio.KEGG import REST
from Bio.KEGG.KGML import KGML_parser
import sys

with open(sys.argv[1]) as f:
    species = f.read().splitlines()

with open(sys.argv[2]) as f:
    list_of_genes = f.read().splitlines()


def listToString(s):  
    
    # initialize an empty string 
    str1 = " " 
    
    # return string   
    return (str1.join(s)) 


def get_sekwencja (id):
   id1='https://www.genome.jp/dbget-bin/www_bget?'+id
   #print(id1)
   r = requests.get(id1)



def  kooo(cccc):
     mydog5=""
     mydog = REST.kegg_find("genes", cccc).read()
     #print(result)
     mydog1 = re.findall('^\S+',mydog)[0]
     #print(mydog1)
     mydog3 = REST.kegg_link("ko", mydog1).read()
     #print("xxx",mydog3)
     if(len(mydog3)<4):
         return(mydog5)
     mydog4=re.findall('ko:\S+',mydog3)[0]
     #print(mydog4)
     mydog5 = REST.kegg_link("genes", mydog4).read()     
     #print(mydog5)
     return(mydog5)



def specie_kegg1(dddd):
    przejscie=re.sub("\nko:\S+",'',dddd);
    mydog6=re.findall('[a-z]+',przejscie)
    #mydog7=re.sub("[\t\s]",'',mydog6)
    return(mydog6)

def gene_kegg(dddd):
    przejscie=re.sub("\nko:\S+",'',dddd);
    mydog6=re.findall('[a-z]+:[a-zA-Z0-9\-_]+',przejscie)    
    #mydog7=re.sub("[\t\s]",'',mydog6)
    return(mydog6)


def dodaj_listy(lst1,lst2):
    df = pd.DataFrame(list(zip(lst1, lst2)), 
               columns =['Name', 'val']) 
    return(df)

def gene_pandas(input):
   return(dodaj_listy(gene_kegg(kooo(input)),specie_kegg1(kooo(input))))


def two_words(input):
    return(re.findall('^\S+[ ]*\S*',str(input))[0])

def longer_name(input):
    p=re.compile('[ ]*\([\S ]+')
    input1=p.sub('',str(input))
    return(input1)

specie=pd.read_pickle('specie')

specie['Organisms.2']=specie['Organisms.1'].apply(two_words)
specie['Organisms.3']=specie['Organisms.1'].apply(longer_name)
#specie.applymap(lambda x: re.findall('^\S+[ ]*\S*',str(x))[0])
specie['Organisms']







#print(gene_pandas('hsa:653604'))
#print(gene_kegg(kooo('hsa:653604')))




def specie_name(name):
    c=0
    c=specie['Organisms'][(specie['Organisms.2']==name)].values[0]
    return(c)



sp1=[]
for sp in species:
    #print(sp,specie_name(sp))
    sp1.append(specie_name(sp))

#print(sp1)



def get_uniprot_id (id):
   id1='https://www.genome.jp/dbget-bin/www_bget?'+id
   #print(id1)
   r = requests.get(id1)
   r1=re.findall('uniprot/[A-Za-z0-9]+',r.text)[0]
   r2=re.sub('uniprot\S','',r1)
   return(r2)

def get_uniprot_lista (lista_id):
   wypluj=''
   for element in lista_id:
       wypluj=wypluj+get_uniprot_id(element)+" "
   return(wypluj)



#tablica_gen_gatunek=gene_pandas('hsa:653604')



#print(tablica_gen_gatunek)


#z=kooo('hsa:653604')
#print(z)
def orthologous_genes(dddd):
   z=kooo(dddd);
   wypluj_lista=[]
   for sp2 in sp1:
     mydog9 = re.findall(sp2+':\S+',z)
     wypluj_lista.append(listToString(mydog9))
   #print(get_uniprot_lista(mydog9))
   #print(get_uniprot_id('lma:LMJF_16_0570K11253'))
   return(wypluj_lista)
  

def append1(input_frame,input_columns,input_values):
          frame_output=pd.DataFrame([input_values],columns=input_values)
          frame_output1=input_frame.append(frame_output)
          return(frame_output1)

species_extend=['gene name']+species
#print(species_extend)
#sp1=[]
def profile_gene1(gene_name):
    #species_extend=['gene name']+species
    pp1=[]
    pp1 += [gene_name]
    profile_extend=pp1+orthologous_genes(gene_name)
    #print(profile_extend)
    frame1=pd.DataFrame([profile_extend],columns=species_extend)
    return(frame1)

def profile_gene_list(list_of_genes):
    #species_extend=['gene name']+species
    #print(species_extend)
    frame0=pd.DataFrame(columns=species_extend)
    for gene in list_of_genes:
         frame0=frame0.append(profile_gene1(str(gene)))
    return(frame0)


wynikowy=profile_gene_list(list_of_genes)


wynikowy=profile_gene_list(list_of_genes)

wynikowy.xls=sys.argv[3]+'.xlsx'
#print(wynikowy.xls)

wynikowy.to_excel(wynikowy.xls,
             sheet_name='Sheet_name_1')
wynikowy.tsv=sys.argv[3]+'.tsv'
wynikowy.to_csv(wynikowy.tsv, sep='\t')

#print(pd.DataFrame([species_extend],columns=species_extend))
#print(profile_gene("YBR248C"))





#print(orthologous_genes('hsa:653604'))
