import os
import pickle
from collections import defaultdict
import grampa_lib.reconcore as RC
import grampa_lib.mul_tree as MT

#############################################################################

#############################################################################

def checkOut(globs):
# Writes output to the checknums file and checks for gene trees over the
# group cap in one or more MUL-trees

	step = "Filtering gene trees over group cap";
	step_start_time = RC.report_step(globs, step, False, "In progress...");
	# Status update

	with open(globs['checknums-outfile'], "w") as checkfile:
		checkfile.write("\t".join(globs['checknums-headers']) + "\n");
		# Open the checknums file and write the headers

		gt_over_cap = defaultdict(list);
		# A dict for the gene trees that are over the cap in one or more mul-trees
		# gt num : [list of mts over cap]

		for mul_num in globs['mul-trees']:
			if mul_num == 0:
				continue;
			# Skip the singly-labeled tree

			cur_mt = globs['mul-trees'][mul_num]
			mt_outline = ["# MT-" + str(mul_num) + ":" + MT.mulPrint(cur_mt[0], cur_mt[2])];
			if not globs['mul-input-flag']:
				mt_outline += ["H1 Node:" + cur_mt[3], "H2 Node:" + cur_mt[4]];
			checkfile.write("\t".join(mt_outline) + "\n");
			# Write some info about the current mul tree as a comment in the checknums file

			groupfilename = os.path.join(globs['pickle-dir'], str(mul_num) + "_groups.pickle");
			cur_groups = pickle.load(open(groupfilename, "rb"));
			# Load the groups from the pickle file

			for gene_num in globs['gt-filtered']:
				if len(globs['gt-filtered'][gene_num]) == 1:
					continue;
				# Skip if the gene tree has already been filtered
				gt_groups, gt_fixed = cur_groups[gene_num][0], cur_groups[gene_num][1];
				num_groups = len(gt_groups);
				num_fixed = len(gt_fixed);
				outline = [mul_num, gene_num, num_groups, num_fixed, 2**num_groups, "N"];
				# Compile the group info for the current gene tree

				if num_groups > globs['cap']:
					gt_over_cap[gene_num].append(mul_num);
					outline[-1] = "Y";
				# Check if the gt is over the group cap for this mt
			
				checkfile.write("\t".join([str(col) for col in outline]) + "\n");
				# Write group info to checknums file for current gt-mt combo
				## End gene tree loop
			checkfile.write("# ----------------------------------\n");
		## End MUL-tree loop

	step_start_time = RC.report_step(globs, step, step_start_time, "Success: " + str(len(gt_over_cap)) + " gts over cap");
	# Status update

	for gt in globs['gt-filtered']:
		if gt in gt_over_cap and len(globs['gt-filtered'][gene_num]) != 1:
			globs['gt-filtered'][gt] = ["# Over group cap in " + str(len(gt_over_cap[gt])) + " MUL-trees"];
			RC.printWrite(globs['logfilename'], globs['log-v'], "# WARNING: Gene tree on line " + str(gt) + " is over the group cap in " + str(len(gt_over_cap[gt])) + " MTs and will be filtered.");
			globs['warnings'] += 1;
	# Change the gts that are over the cap to a filter string in the list of gts

	return globs, len(gt_over_cap);
	# This block handles output to the checknums file and filters any trees over the group cap (-c).
	# The call of the important collapseGroups function that groups polyploid clades in the gene trees to speed up the reconciliations to MUL-trees.

#############################################################################

def filterOut(globs, num_filtered_gt):

	gene_trees_final = {};
	if num_filtered_gt != 0:
		RC.printWrite(globs['logfilename'], globs['log-v'], "# " + RC.getDateTime() + "  INFO: Filtered gene trees: ", str(num_filtered_gt), 40);
		RC.printWrite(globs['logfilename'], globs['log-v'], "# " + RC.getDateTime() + "  INFO: Final gene tree count for reconciliations: ", str(len(globs['gt-filtered']) - num_filtered_gt), 40);
		# Some info

		step = "Writing filtered gene trees to file";
		step_start_time = RC.report_step(globs, step, False, "In progress...");
		# Status update
		
		gt_written = 0;
		with open(globs['gt-filtered-file'], "w") as filterfile:
			for gene_num in globs['gt-filtered']:
				filterfile.write(globs['gt-filtered'][gene_num][0] + "\n");
				if len(globs['gt-filtered'][gene_num]) != 1:
					gene_trees_final[gene_num] = globs['gt-filtered'][gene_num];
					gt_written += 1;
		# Write gene trees to a file if they aren't filtered

		step_start_time = RC.report_step(globs, step, step_start_time, "Success: " + str(gt_written) + " gene trees written");
		# Status update
	## Filter block

	else:
		RC.printWrite(globs['logfilename'], globs['log-v'], "# " + RC.getDateTime() + "  INFO: No trees filtered! Using your original set.");
		gene_trees_final = globs['gt-filtered'];
	## If no gene trees are filtered, update here

	globs['gt-filtered'] = gene_trees_final;
	return globs;

#######################################################['mul.tree', 'h1.node', 'h2.node', 'labeled.tree', 'score']######################

def mainOut(globs, all_scores, lowest_maps, multiple_maps):

	if globs['mul-input-flag']:
		globs['output-headers'] = [globs['output-headers'][0]] + globs['output-headers'][3:];
	# Adjust the headers if the input was a MUL-tree

	with open(globs['output-file'], "w") as outfile, open(globs['dup-count-outfile'], "w") as dupout:
		outfile.write("\t".join(globs['output-headers']) + "\n");
		dupout.write("\t".join(globs['dup-count-headers']) + "\n");
		# Write the headers

		count = 0;
		for mul_num, mul_score in all_scores:
			mul_tree = globs['mul-trees'][mul_num];

			mt, hybrid_clade, hybrid_node, copy_node = mul_tree[0], mul_tree[2], mul_tree[3], mul_tree[4];
			# Unpack the current mul tree

			if mul_num == 0:
				outline = [mul_num, "NA", "NA", mul_score, globs['parsed-st-str']];
			# For the input singly-labeled tree, no h nodes are output
			else:
				if globs['mul-input-flag']:
					outline = [mul_num, mul_score, MT.mulPrint(mt, hybrid_clade)];
				else:
					outline = [mul_num, hybrid_node, copy_node, mul_score, MT.mulPrint(mt, hybrid_clade)];
			# Output for the mul trees

			outfile.write("\t".join([str(col) for col in outline]) + "\n");
			# Write to the main output

			if count < 6:
				min_dict = mul_tree[1];
				# Unpack the dictionary of the mul tree

				main_dups = { node : 0 for node in min_dict };
				# Dict to count dups per node

				for m in lowest_maps[count]:
					maps, dups = lowest_maps[count][m][0][3], lowest_maps[count][m][0][4];
					# Unpack the maps for the current mul tree

					for gt_node in dups:
						if dups[gt_node] != 0:
							main_dups[maps[gt_node][0]] += 1;	
					# Count dups per species tree node for every gene tree
				## End map loop

				for node in main_dups:
					outline = [str(mul_num), node, str(main_dups[node])];
					dupout.write("\t".join(outline) + "\n");			
				# Write out the counts per node
			## For the lowest scoring trees, also output counts of duplications per node

			count += 1;
		## End MUL-tree loop

	# if multiple_maps != 0:
	# 	num_str = " trees have ";
	# 	if multiple_maps == 1:
	# 		num_str = " tree has ";
	# 	multiple_outline = "# MSG: " + str(multiple_maps) + " gene" + num_str + "multiple maps to the species tree with equal scores. Only one of these maps is (randomly) chosen in the final duplication counts. See detailed output file for more info."
	# 	RC.printWrite(globs.dupcountfilename, globs.main_v, "# ----------------------------------------");
	# 	RC.printWrite(globs.dupcountfilename, globs.main_v, multiple_outline)

	## Close the files

	return;

		# # RC.printWrite(globs.outfilename, globs.main_v, "# ----------------------------------------");
		# # if globs.spec_type == 's':
		# # 	if min_num != 0:
		# # 		RC.printWrite(globs.outfilename, globs.v, "# The MUL-tree with the minimum parsimony score is MT-" + str(min_num) + ":\t" \
		# # 			+ MT.mulPrint(min_tree[0], min_tree[2]));
		# # 	else:
		# # 		RC.printWrite(globs.outfilename, globs.v, "# The tree with the minimum parsimony score is the singly-labled tree (ST):\t" + min_tree[0]);
		# # 	RC.printWrite(globs.outfilename, globs.v, "# Score = " + str(min_score));
		# # 	RC.printWrite(globs.outfilename, globs.main_v, "# ----------------------------------------");

#############################################################################

def detOut(globs, min_maps):
	
	if not globs['maps-opt']:
		globs['detailed-headers'] = globs['detailed-headers'][:-1];
	else:
		import grampa_lib.gene_tree as GT
	# Remove the maps header if --maps isn't set
	# Otherwise, import the gene_tree library to parse the maps

	multiple_maps = 0;

	with open(globs['detailed-outfile'], "w") as detfile:
		detfile.write("\t".join(globs['detailed-headers']) + "\n");
		# Write the headers for the detailed output file

		for gene_num, cur_maps in min_maps.items():

			mul_map_string = "";
			if len(cur_maps) != 1:
				multiple_maps += 1;

				mulmap_outline = "# GT-" + str(gene_num) + " to MT-" + str(globs['min-num']) + "\t" + str(len(cur_maps)) + " maps found!"
				detfile.write(mulmap_outline + "\n")

			for cur_map in cur_maps:
				outline = [globs['min-num'], gene_num, cur_map[1], cur_map[2], cur_map[0]];

				if globs['maps-opt']:
					outline.append(GT.detailedOut(globs['gt-filtered'][gene_num][0], globs['gt-filtered'][gene_num][1], cur_map[3], cur_map[4], cur_map[5]));

				detfile.write("\t".join([str(col) for col in outline]) +"\n");
				#RC.printWrite(globs.detoutfilename, globs.v, outline);
	
	return multiple_maps;
	# Output to the detailed file.

#############################################################################