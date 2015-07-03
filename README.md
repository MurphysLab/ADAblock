ADAblock for ImageJ
--------------------
***Automated Defect Analysis of Block Copolymers***


*Jeffrey N. Muprhy*
[@MurphysLab](https://twitter.com/MurphysLab)

---------------------

Code written for the analysis of line patterns produced from the domains of block copolymer thin films.
Sample SEM images can be obtained from [https://hdl.handle.net/10402/era.41438](https://hdl.handle.net/10402/era.41438).

The code is made available for use and modification under the "[MIT License](https://opensource.org/licenses/MIT)".

If you have any questions, or wish to contribute to the code, please contact the author.

---------------------

**Scripts:**

* **ADAblock.ijm** — The algorithm is designed to perform analysis of images of BCP thin films or surfaces structured by BCP thin films. It is run as a macro within ImageJ. 
* **Data_Amalgamation_Script.py** — Collects data from each of the output files produced by ADAblock.ijm, in the directories produced by ADAblock, and collects them into a single CSV file.
* **Data_Sort_and_Plot.py** — The data sort and plot script takes data from the CSV file produced by Data_Amalgamation_Script.py and can be used to re-sort and plot the data.

---------------------

Download ImageJ: https://rsbweb.nih.gov/ij

Download Python: https://www.python.org

Open a 7-zip archive: https://www.7-zip.org