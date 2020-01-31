 # Protein-Protein Interfaces and Binding Site Visualization

 In many protein interactions, structures which consist of 2 or more interacting proteins , called protein-protein complexes are observed. This project involves the determination of the amino acids which are present at the interface of two interacting protein chains.

 The interface of interacting protein chains consists of a set of amino acids belonging to the two protein chains within a certain threshold distance. The script calculates the Euclidean distance between the amino acids and determines the interface atoms.

 **Input**:

      The following need to specified as command-line arguments while running the script:

          * path for .pdf format file of a particular protein (-f)
          * threshold value for calculating and specifying the Euclidean distance between the interface atoms (-t)
          * Name of the first chain (Should be a single capital alphabet) (-c1)
          * Name of the second chain (Should be a single capital alphabet) (-c2)

**Code Execution:** ./protein_binding_site.py -f <pdb_file> -t <threshold value> -c1 <chain1> -c2 <chain2>

**Example Output for:**
