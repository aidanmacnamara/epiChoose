# epiChoose
Data-driven selection of cell line models and epigenetic comparison across cell types

* 12/07/2017

	* Current matrix lists

		* all_data - Blueprint + Project 1 data + some lung cell line outliers, size = 50 x 18068, generated from vignettes/project_1.R, TSSs of protein-coding genes
	
		* mask_data - Blueprint + Project 1/3 + ENCODE, size = 112 x 444967, generated from test/test_masking.R, all regulatory regions (defined by Ensembl multi-cell build) 

		* mask_data (Aidan, local) - As above with Project 2 data added, size = 137 x 444967
	
	Current status
	
		* Graeme (Tessella) - working with mask_data using a combination of PCA and PLS (?) to attempt to uncover an immortalisation signature
		
		* Nam - working with all_data using hierarchical clustering on combined data types to model the haematopoietic lineage hieraerchy
		
		* David - Using gene/datatype clustering to look for differential patterns between cell groups
		