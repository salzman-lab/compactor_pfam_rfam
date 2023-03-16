# compactor_pfam_rfam
A python script, 2 bash scripts, and a Sherlock sbatch to run Rfam and Pfam on Adam Gudys's compactor tabular output, then to report the best Rfam and Pfam hits for each compactor in a new table.



The best Rfam and Pfam hits are reported by e-value in Rfam and full sequence e-value in Pfam. 
A compactor is selected for Pfam or Rfam alignment if exact_support > 0. 

The user needs only to specify the compactor output file and a run name, then to type "sbatch main.sh" to run this code in the Sherlock computing environment. 
