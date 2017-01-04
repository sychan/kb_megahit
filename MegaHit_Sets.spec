/*
A KBase module to wrap the MEGAHIT package.
*/

module MegaHit_Sets {

	/* run_megahit()

	    Run MEGAHIT.  Most parameters here are just passed forward to MEGAHIT

	    run_megahit() is responsible for accepting input params from Narrative, 
	        calling exec_megahit(), and generating report.  
		It mediates communication with the Narrative

	    workspace_name - the name of the workspace for input/output
	    input_reads_ref - the ref to the PE read library or ReadsSet
	                       (SE library support in the future)
	    output_contig_set_name - the base name of the output contigset or AssemblySet

	    combined_assembly_flag - if input is a ReadsSet, indicate combined Assembly

	    megahit_parameter_preset - override a group of parameters; possible values:
	        
		meta            '--min-count 2 --k-list 21,41,61,81,99'
                                (generic metagenomes, default)
		meta-sensitive  '--min-count 2 --k-list 21,31,41,51,61,71,81,91,99'
                                 (more sensitive but slower)
	        meta-large      '--min-count 2 --k-list 27,37,47,57,67,77,87'
                                (large & complex metagenomes, like soil)
		bulk            '--min-count 3 --k-list 31,51,71,91,99 --no-mercy'
                                (experimental, standard bulk sequencing with >= 30x depth)
		single-cell     '--min-count 3 --k-list 21,33,55,77,99,121 --merge_level 20,0.96'
                                (experimental, single cell data)

	     min_count - minimum multiplicity for filtering (k_min+1)-mers, default 2

    	     min_k - minimum kmer size (<= 127), must be odd number, default 21
    	     max_k - maximum kmer size (<= 127), must be odd number, default 99
             k_step - increment of kmer size of each iteration (<= 28), must be even number, default 10
	     
             k_list - list of kmer size (all must be odd, in the range 15-127, increment <= 28);
                 overrides '--k-min', '--k-max', and '--k-step'

	     min_contig_length - minimum length of contigs to output, default 200
        */

	/* Kmer Params
	**     @optional min_count
	**     @optional k_min
	**     @optional k_max
	**     @optional k_step
	**     @optional k_list
	*/
	typedef structure {
		int min_count;
		int k_min;
		int k_max;
		int k_step;
		list <int> k_list;
	} Kmer_Params;	    

	/* run_megahit()
	**
	**     @optional megahit_parameter_preset
	**     @optional min_contig_len
	*/
	typedef structure {
		string workspace_name;
		string input_reads_ref;
		string output_contigset_name;
		int combined_assembly_flag;  /* 1=True, 0=False, def: 1 */
		string megahit_parameter_preset;

		int min_contig_len;
		Kmer_Params kmer_params;
	} MegaHitParams;

	typedef structure {
		string report_name;
		string report_ref;
	} MegaHitOutput;

	funcdef run_megahit(MegaHitParams params) returns (MegaHitOutput output)
		authentication required;


	/* exec_megahit()

	    Actual execution of MEGAHIT!

	    Accepts ReadsSet or a ReadsLibrary as Input

	    Creates Assembly object(s) as output.
	    Will eventually also create AssemblySet object if input is a ReadsSet and not running a combined assembly
	
	    Other vars same as run_megahit()
	*/
	typedef structure {
		string workspace_name;
		string input_reads_ref;
		string output_contigset_name;
		int combined_assembly_flag;  /* 1=True, 0=False, def: 1 */
		string megahit_parameter_preset;

		int min_count;
		int k_min;
		int k_max;
		int k_step;
		list <int> k_list;
		int min_contig_len;
	} ExecMegaHitParams;

	typedef structure {
	        string       report_text;
		list<string> output_contigset_ref;
	} ExecMegaHitOutput;

	funcdef exec_megahit(ExecMegaHitParams params) returns (ExecMegaHitOutput output)
		authentication required;
};