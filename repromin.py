def combinatorial_analysis(network, proteomic_data, TF_list, n_combinations, export_file_name, units="fg"):
	''''''
		
	import pandas as pd
	import numpy as np
	import math
	import itertools as it
	import time
	
	##### LOAD DATA ####

	a=pd.read_csv(network)
	a.columns = ["TF","Target","Regulation"]
	b=pd.read_csv(proteomic_data)
	b.columns = ["Target","PB"]
	c=a.merge(b,on="Target", how="left").drop_duplicates().fillna(0)

	#### VARIABLES ####
	
	start_time = time.time()
	n = 0
	total_proteome = b["PB"].values.sum()
	top_proteomic_balance = 0
	results=pd.DataFrame(columns = ["Combination", "# of Silenced", "# of Induced", "PB " + units, " Silenced Targets", "Induced Targets"])
	all_combinations = list(it.combinations(TF_list,n_combinations))
	
	print(len(all_combinations), " combinations to test")
	
	##### COMBINATORIAL ######
	
	for combination in all_combinations:
		
		n = n + 1
		all_pos_tar, all_neg_tar = [], []
		i, j = [], []

		for tf in combination:

			pos_tar=list(c.loc[(c["TF"] == tf) & (c["Regulation"] == "+")].Target.values)
			neg_tar=list(c.loc[(c["TF"] == tf) & (c["Regulation"] == "-")].Target.values)

			all_pos_tar=all_pos_tar+pos_tar
			all_neg_tar=all_neg_tar+neg_tar

		all_pos_tar=pd.Series(all_pos_tar)
		all_pos_tar=all_pos_tar.value_counts()

		all_neg_tar=pd.Series(all_neg_tar)
		all_neg_tar=all_neg_tar.value_counts()


		for target, pos_int in zip (all_pos_tar.index, all_pos_tar.values): 
			pos_regulators = c.loc[ (c["Target"] == target) & (c["Regulation"] == "+") ].Target.count()-pos_int
			neg_regulators = c.loc[ (c["Target"] == target) & (c["Regulation"] == "-") ].Target.count()
			
			if (pos_regulators == 0) & (neg_regulators == 0):
				i.append(target)

		silenced = pd.DataFrame()
		silenced["Target"] = i
		silenced = silenced.merge(c.iloc[:,[1,3]], on="Target").drop_duplicates()
		silenced_pb = silenced["PB"].sum()

		if math.isnan(silenced_pb) == True:
				silenced_pb = 0

		for target, neg_int in zip (all_neg_tar.index, all_neg_tar.values):
			pos_regulators=c.loc[ (c["Target"] == target) & (c["Regulation"] == "+") ].Target.count()
			neg_regulators=c.loc[ (c["Target"] == target) & (c["Regulation"] == "-") ].Target.count()-neg_int
			
			if (neg_regulators == 0) & (pos_regulators == 0):
				j.append(target)

		induced = pd.DataFrame()
		induced["Target"] = j
		induced = induced.merge(c.iloc[:,[1,3]], on="Target").drop_duplicates()
		induced_pb = induced["PB"].sum()

		if math.isnan(induced_pb) == True:
				induced_pb = 0

		proteomic_balance = silenced_pb - induced_pb
		results.loc[n]=[combination,len(i),len(j),proteomic_balance,i,j]
		
		if proteomic_balance > top_proteomic_balance:
			
			top_combination = combination
			top_silenced = len(i)
			top_induced = len(j)
			top_proteomic_balance = proteomic_balance
			top_silenced_genes = i
			top_induced_genes = j 
		
	results.to_csv(export_file_name)
	elapsed_time = time.time() - start_time
	print(" The Best Combination is:")
	print(top_combination, " - ", top_proteomic_balance, " ", units, " liberated. ", "Corresponding to: ", top_proteomic_balance/total_proteome*100, "% of total proteome.")    
	print("Elapsed Time: ", elapsed_time/60, " minutes")
	
def TRN_analysis(essential_genes, network, proteomic_data, essential_TF_export_filename, nonessential_TF_export_filename, units="fg", refined=False):
	''' '''
	import pandas as pd
	import numpy as np
	import math

	#### LOAD DATA ####

	a=pd.read_csv(network)
	a.columns = ["TF","Target","Regulation"]
	b=pd.read_csv(proteomic_data)
	b.columns = ["Target","PB"]
	d=pd.read_csv(essential_genes)
	d.columns = ["Target"]
	d["Essential"] = ["Yes"]*len(d["Target"])
	c=a.merge(b,on="Target", how="left").drop_duplicates().fillna(0)
	c=c.merge(d,on="Target", how="left").drop_duplicates().fillna("No")

	#### TF ESSENTIALITY CLASIFICATION ###

	n=0
	tfs=np.unique(c["TF"])
	TF_targets=pd.DataFrame(columns=["TF","Positive Essential"])
	essential_TFs=pd.DataFrame()
	
	for tf in tfs:
		
		n=n+1
		pos_tar= c.loc[(c["TF"] == tf) & (c["Regulation"] == "+") & (c["Essential"] == "Yes")].Target.count() 
		
		TF_targets.loc[n] = [tf,pos_tar]
					 
	ne_tfs=TF_targets[TF_targets["Positive Essential"] == 0].TF.values
	es_tfs=TF_targets[TF_targets["Positive Essential"] > 0].TF.values
	essential_TFs["Essential TF"]=es_tfs
	essential_TFs.to_csv(essential_TF_export_filename)

	#### CALCULATE PROTEOMIC BALANCE OF NON-ESSENTIAL TF ###

	real_pb=pd.DataFrame()
	x,y,z,zz=[],[],[],[]

	for tf in ne_tfs:
		
		pos_tar=list(c[(c["TF"] == tf) & (c["Regulation"] == "+")].Target.values)
		neg_tar=list(c[(c["TF"] == tf) & (c["Regulation"] == "-")].Target.values)

		i,j,k,l,m,n,=[],[],[],[],[],[]
		pos_ex_target=pd.DataFrame()
		neg_ex_target=pd.DataFrame()

		for target in pos_tar:
			pos_regulators= c.loc[ (c["Target"] == target) & (c["Regulation"] == "+") ].Target.count()
			neg_regulators= c.loc[ (c["Target"] == target) & (c["Regulation"] == "-") ].Target.count()
				  
			if (pos_regulators == 1) & (neg_regulators == 0):
				i.append(target)
				j.append(pos_regulators)
				k.append(neg_regulators)
					   
		pos_ex_target["Target"]=i
		pos_ex_target["Positive Regulators"]=j
		pos_ex_target["Negative Regulatos"]=k
		pos_ex_target=pos_ex_target.merge(c.iloc[:,[1,3]],on="Target")
		pos_ex_target=pos_ex_target.drop_duplicates()
		
		pos_pb=pos_ex_target["PB"].sum()
		pos_count=pos_ex_target["PB"].count()
		
		if math.isnan(pos_pb)== True:
			pos_pb=0
		if math.isnan(pos_count)== True:
			pos_count=0
		
		x.append(pos_pb)
		y.append(pos_count)

		
		for target in neg_tar:
			pos_regulators=c.loc[ (c["Target"] == target) & (c["Regulation"] == "+") ].Target.count()
			neg_regulators=c.loc[ (c["Target"] == target) & (c["Regulation"] == "-") ].Target.count()
		   
			if (neg_regulators == 1) & (pos_regulators == 0):
				l.append(target)
				m.append(pos_regulators)
				n.append(neg_regulators)

		neg_ex_target["Target"]=l
		neg_ex_target["Positive Regulators"]=m
		neg_ex_target["Negative Regulatos"]=n
		neg_ex_target=neg_ex_target.merge(c.iloc[:,[1,3]],on="Target")
		neg_ex_target=neg_ex_target.drop_duplicates()
		
		neg_pb=neg_ex_target["PB"].sum()
		neg_count=neg_ex_target["PB"].count()
		
		if math.isnan(neg_pb) == True:
			neg_pb=0
		if math.isnan(neg_count) == True:
			neg_count=0
		
		z.append(neg_pb)
		zz.append(neg_count)

	real_pb["TF"]=ne_tfs
	real_pb["Positive Count"]=y
	real_pb["Negative Count"]=zz
	real_pb["Positive PB"]=x
	real_pb["Negative PB"]=z
	real_pb["Real PB"]=real_pb["Positive PB"]-real_pb["Negative PB"]
	real_pb.sort_values("Real PB", ascending = 'False')
	real_pb.to_csv(nonessential_TF_export_filename)

	#### DEFINE CANDIDATE TF ###
	
	candidate_TF = real_pb[(real_pb["Positive Count"] > 0)].TF.values
	positive_candidate_TF = real_pb[(real_pb["Positive Count"] > 0) & (real_pb["Real PB"] > 0)].TF.values
	
	print("Essential TF list and Non-essential TF Proteomic Balance was exported.\n")
	
	if not refined:
		print("Candidate TFs:")
		print(candidate_TF)
		return candidate_TF
	
	if refined:
		print("Candidate TFs with positive PB:")
		print(positive_candidate_TF)
		return positive_candidate_TF

def search_TF_combination(network, proteomic_data, TF_combination, units="fg"):
	''''''
	import pandas as pd
	import numpy as np
	import math
		
	##### LOAD DATA ####

	a=pd.read_csv(network)
	a.columns = ["TF","Target","Regulation"]
	b=pd.read_csv(proteomic_data)
	b.columns = ["Target","PB"]
	c=a.merge(b,on="Target", how="left").drop_duplicates().fillna(0)

	#### VARIABLES ####
	
	all_pos_tar, all_neg_tar = [], []
	i, j = [], []
	n = 0
	total_proteome = b["PB"].values.sum()

	for tf in TF_combination:
		
		if tf in c["TF"].values:

			pos_tar=list(c.loc[(c["TF"] == tf) & (c["Regulation"] == "+")].Target.values)
			neg_tar=list(c.loc[(c["TF"] == tf) & (c["Regulation"] == "-")].Target.values)

			all_pos_tar=all_pos_tar+pos_tar
			all_neg_tar=all_neg_tar+neg_tar
			n=n+1			

		else:
			
			print("Warning ",tf," not found") 
			
	all_pos_tar=pd.Series(all_pos_tar)
	all_pos_tar=all_pos_tar.value_counts()

	all_neg_tar=pd.Series(all_neg_tar)
	all_neg_tar=all_neg_tar.value_counts()


	for target, pos_int in zip (all_pos_tar.index, all_pos_tar.values): 
		pos_regulators = c.loc[ (c["Target"] == target) & (c["Regulation"] == "+") ].Target.count()-pos_int
		neg_regulators = c.loc[ (c["Target"] == target) & (c["Regulation"] == "-") ].Target.count()
			
		if (pos_regulators == 0) & (neg_regulators == 0):
			i.append(target)

	silenced = pd.DataFrame()
	silenced["Target"] = i
	silenced = silenced.merge(c.iloc[:,[1,3]], on="Target").drop_duplicates()
	silenced_pb = silenced["PB"].sum()

	if math.isnan(silenced_pb) == True:
			silenced_pb = 0

	for target, neg_int in zip (all_neg_tar.index, all_neg_tar.values):
		pos_regulators=c.loc[ (c["Target"] == target) & (c["Regulation"] == "+") ].Target.count()
		neg_regulators=c.loc[ (c["Target"] == target) & (c["Regulation"] == "-") ].Target.count()-neg_int
			
		if (neg_regulators == 0) & (pos_regulators == 0):
			j.append(target)

	induced = pd.DataFrame()
	induced["Target"] = j
	induced = induced.merge(c.iloc[:,[1,3]], on="Target").drop_duplicates()
	induced_pb = induced["PB"].sum()

	if math.isnan(induced_pb) == True:
				induced_pb = 0

	proteomic_balance = silenced_pb - induced_pb	
	
	print("Result: ",TF_combination, " - ", proteomic_balance, " ", units, " liberated. ", "Corresponding to: ",proteomic_balance/total_proteome*100,"% of total proteome.")
	print("Number of TFs counted: ", n, "\n")   
	