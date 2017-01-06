# iMM904
Implementing MOMA algorithm into PSAMM. Validating results with the implementation results from <a href="http://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-3-37">this paper</a>.
# Important Scripts
<b>moma.py</b> : 4 Implementations of the MOMA algorithm (LP2, QLP2, LP3, QLP3)<br>
<b>moma_test.py</b> : The script that deletes the genes and reports the percent viability compared to the wild type<br>
<b>gene_deletion.py</b> : The script sets up the environment for each of the media experiments. Ultimately calls the moma_test.py script<br>

# Important Functions
<b>moma_test.gene_deletion(model, mm, solver, algorithm)</b>: algorithm is a parameter to select which MOMA algorithm to use<br>
