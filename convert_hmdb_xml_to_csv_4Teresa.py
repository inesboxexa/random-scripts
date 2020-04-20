print ('Importing smorgasbord...')
#from pathlib import Path 
import lxml.etree as ET
import pandas as pd
from tqdm import tqdm

tag='{http://www.hmdb.ca}hmdb'

#tree = ET.parse('example.xml')
context = ET.iterparse('hmdb_metabolites_v2.xml', events=("end",))
_, root = next(context)  # Grab the root element.

data={}
print('Running...')
for event, element in  tqdm(context):
	if element.tag=='{http://www.hmdb.ca}metabolite':		 
		hmdb_id={}
		name={}
		monisotopic_molecular_weight={}
		kingdom={}
		super_class={}
		classs={}
		sub_class={}

		for children in element.getchildren():
			if children.tag=="{http://www.hmdb.ca}accession":
				hmdb_id=children.text

			if children.tag=="{http://www.hmdb.ca}name":
				name["name"]=children.text

			if children.tag=="{http://www.hmdb.ca}monisotopic_molecular_weight":
				monisotopic_molecular_weight["monisotopic_molecular_weight"]=children.text
				
			if children.tag=="{http://www.hmdb.ca}taxonomy":
				for taxonomy_children in children.getchildren():			
					if taxonomy_children.tag=="{http://www.hmdb.ca}kingdom":						
						kingdom["kingdom"]=taxonomy_children.text
			
				for taxonomy_children in children.getchildren():			
					if taxonomy_children.tag=="{http://www.hmdb.ca}super_class":						
						super_class["super_class"]=taxonomy_children.text
					
				for taxonomy_children in children.getchildren():			
					if taxonomy_children.tag=="{http://www.hmdb.ca}class":						
						classs["class"]=taxonomy_children.text
				for taxonomy_children in children.getchildren():			
					if taxonomy_children.tag=="{http://www.hmdb.ca}sub_class":						
						sub_class["sub_class"]=taxonomy_children.text		

		data[hmdb_id]=name,monisotopic_molecular_weight,kingdom,super_class,classs,sub_class
	root.clear()	                                

print("Converting to Dataframe")
hmdb_id = pd.DataFrame([])
name = pd.DataFrame([])
monisotopic_molecular_weight = pd.DataFrame([])
kingdom = pd.DataFrame([])
super_class = pd.DataFrame([])
classs = pd.DataFrame([])
sub_class = pd.DataFrame([])

for key,value in tqdm(data.items()):
    for items in value:
        for key1,value1 in items.items():
            if key1=="name":
                name=name.append(pd.DataFrame({"name":[value1]},index=[key])) 
            elif key1=="monisotopic_molecular_weight":
                 monisotopic_molecular_weight  =monisotopic_molecular_weight  .append(pd.DataFrame({"monisotopic_molecular_weight":[value1]},index=[key]))
            elif key1=="kingdom":                
                 kingdom=kingdom.append(pd.DataFrame({"kingdom":[value1]},index=[key]))
            elif key1=="super_class": 
                 super_class =super_class.append(pd.DataFrame({"super_class":[value1]},index=[key]))
            elif key1=="class": 
                 classs =classs.append(pd.DataFrame({"class":[value1]},index=[key]))
            elif key1=="sub_class": 
                 sub_class =sub_class.append(pd.DataFrame({"sub_class":[value1]},index=[key]))

result = pd.concat([hmdb_id,name,monisotopic_molecular_weight,kingdom,super_class,classs,sub_class ],axis=1)
result.to_csv("hmdb_metabolites_teresa2.csv")

print('End without errors')