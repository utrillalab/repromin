import pandas as pd
import numpy as np
import math
from itertools import combinations
from time import time
import os.path

def combinatorial_analysis(network_data, proteomic_data, TF_list, n_combinations, export_file_name, full_export=True, units="fg", verbose=False):
    
    print("Analysis in Progress. Preparing Data.")
    
    # Global Data
    start_time = time()
    network = pd.read_csv(network_data)
    proteomics = pd.read_csv(proteomic_data)
    total_proteome = proteomics.iloc[:, 1].values.sum()
    positive_targets, negative_targets, targets_positive_degree, targets_negative_degree, targets_proteomics = {}, {}, {}, {}, {}
    log_path = export_file_name+"_Elapsed_Times.csv"
    
    # Create Log file if it does not already exist
    if not os.path.isfile(log_path):
        log = pd.DataFrame(columns = ("#KO", "Elapsed_Time(seg)", "#Combinations", "Seed_Size"))
        log.to_csv(log_path)

    # Index TFs Targets
    for TF in network.iloc[:, 0].unique():  
        positive_targets[TF] = network[(network.iloc[:, 0] == TF) & (network.iloc[:, 2] == "+")].iloc[:,1].values.tolist()
        negative_targets[TF] = network[(network.iloc[:, 0] == TF) & (network.iloc[:, 2] == "-")].iloc[:,1].values.tolist()

    # Index Targets In-degree
    for target in network.iloc[:, 1].unique():
        targets_positive_degree[target] = len(network[(network.iloc[:, 1] == target) & (network.iloc[:, 2] == "+")])
        targets_negative_degree[target] = len(network[(network.iloc[:, 1] == target) & (network.iloc[:, 2] == "-")])

    # Index Targets Proteomics
    for gene in proteomics.iloc[:, 0]:
        targets_proteomics[gene] = proteomics[proteomics.iloc[:, 0] == gene].iloc[0,1]
    
    elapsed_time = time() - start_time
    print("Elapsed Time Preparing Data:", elapsed_time, "seconds\n")

    # Start Combinatorial
    for n in n_combinations:
        
        print("Solving", n, "KO")
        start_time = time()
        number_of_combinations = int(math.factorial(len(TF_list))/math.factorial(n)/math.factorial(len(TF_list)-n))
        top_proteomic_balance = 0
        results = {}
        counter, progress = 0, 0
        progress_bar = [int(number_of_combinations*n) for n in np.arange(0, 1.1, .1)]
        
        print(number_of_combinations, " combinations to test\n")
        
        # Calculate All N-Combinations
        for combination in combinations(TF_list, n):
            combination_pos_targets, combination_neg_targets, silenced, induced = [], [], [], []
            silenced_pb, induced_pb = 0, 0
            set_flag = True

            # Find Combination Targets
            for TF in combination:
                combination_pos_targets = combination_pos_targets + positive_targets[TF]
                combination_neg_targets = combination_neg_targets + negative_targets[TF]

            # Find Combination Silenced Genes 
            for pos_target in combination_pos_targets:
                if (combination_pos_targets.count(pos_target) == targets_positive_degree[pos_target]) & (targets_negative_degree[pos_target] == 0):
                    silenced.append(pos_target)
            
            # Find Combination Induced Genes 
            for neg_target in combination_neg_targets:
                if (combination_neg_targets.count(neg_target) == targets_negative_degree[neg_target]) & (targets_positive_degree[neg_target] == 0):
                    induced.append(neg_target)
            
            # Calculate Combination Proteomic Balance
            for gene in silenced:
                if gene in targets_proteomics.keys():
                    silenced_pb = silenced_pb + targets_proteomics[gene]

            for gene in induced:
                if gene in targets_proteomics.keys():
                    induced_pb = induced_pb + targets_proteomics[gene]

            proteomic_balance = silenced_pb - induced_pb
            percentaje = proteomic_balance/total_proteome*100

            # Check if Combination is Top or Low
            if proteomic_balance >= top_proteomic_balance:
                top_combination = combination
                top_silenced = len(silenced)
                top_induced = len(induced)
                top_proteomic_balance = proteomic_balance
                top_silenced_genes = silenced
                top_induced_genes = induced
                top_percentaje = top_proteomic_balance/total_proteome*100

            if set_flag:
                low_proteomic_balance = top_proteomic_balance
                low_combination = top_combination
                low_silenced = top_silenced
                low_induced = top_induced
                low_silenced_genes = top_silenced_genes
                low_induced_genes = top_induced_genes
                low_percentaje = top_percentaje
                set_flag = False

            if proteomic_balance < low_proteomic_balance:
                low_combination = combination
                low_silenced = len(silenced)
                low_induced = len(induced)
                low_proteomic_balance = proteomic_balance
                low_silenced_genes = silenced
                low_induced_genes = induced
                low_percentaje = low_proteomic_balance/total_proteome*100
            
            # Save all results
            if full_export:
                results[len(results)] = [combination, len(silenced), len(induced), proteomic_balance, percentaje, silenced, induced]
            
            # Track Progress
            if verbose:
                counter = counter + 1
                
                if counter in progress_bar:
                    progress = progress + 10
                    elapsed_time = time() - start_time
                    print("Progress:", progress,"%")
                    print("Elapsed Time", elapsed_time//60,"min")
                
        # Save only "Best" and "Worst" Combination
        print("\n")
        if not full_export:
            results["Best"] = [top_combination, top_silenced, top_induced, top_proteomic_balance, top_percentaje, top_silenced_genes, top_induced_genes]
            results["Worst"] = [low_combination, low_silenced, low_induced, low_proteomic_balance, low_percentaje, low_silenced_genes, low_induced_genes]

        elapsed_time = time() - start_time

        print("The Best Combination is:", top_combination)
        print("Total of:", top_proteomic_balance, units, "liberated.\nCorresponding to ", top_proteomic_balance/total_proteome*100, "% of total proteome.\n")
        
        if verbose:
            print("The Worst Combination is:", low_combination)
            print("Total of:", low_proteomic_balance, units, "liberated.\nCorresponding to ", low_proteomic_balance/total_proteome*100, "% of total proteome.\n")   

        print("Elapsed Time in Combinatorial Analysis:", elapsed_time, "seconds")    

        # Generate Export File
        export = pd.DataFrame(results).T
        export.columns = ["Combination", "# of Silenced", "# of Induced", "PB " + units, "%Proteome", "Silenced Targets", "Induced Targets"]
         
        log = pd.read_csv(log_path, index_col = 0)
        log.loc[len(log)] = [n, elapsed_time, number_of_combinations, len(TF_list)]
        
        if full_export:
            export.to_csv(export_file_name+"_"+str(n)+"KO.csv")
            log.to_csv(log_path)
            print("Full Results Exported Succesfully!")
        else:
            export.to_csv(export_file_name+"_"+str(n)+"KO_Partial.csv")
            log.to_csv(log_path)
            print("Partial Results Exported Succesfully!")
            
        print("\n")

def TRN_analysis(essential_genes, network, proteomic_data, export_file_name1, export_file_name2, units="fg", refined=False):

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
	essential_TFs.to_csv(export_file_name1)

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
	real_pb.to_csv(export_file_name2)

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
	print("Silenced Genes:", len(i), i)
	print("Induced Genes:", len(j), j)
	print("Number of TFs counted: ", n, "\n")
    
def get_gene_names(bnumbers, gene_database_path):
    
    df = pd.read_csv(gene_database_path)
    a = []
    x = pd.DataFrame()
    
    for b in bnumbers:
        
        if len(df[df["b"] == b]["n"]) == 1:
            c = df[df["b"] == b]["n"].values[0]
            a.append(c)
        
        else:
            print("Gene: ", b, "Not Found")
    
    x["Gene"] = a
    
    return x