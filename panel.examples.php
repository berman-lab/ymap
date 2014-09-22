<?php
	session_start();
	if (isset($_SESSION['logged_on'])) { $user = $_SESSION['user']; } else { $user = 'default'; }
?>
<style type="text/css">
	html * {
		font-family: arial !important;
	}
</style>
<font size='3'>Datasets used in associated paper.</font><br><br>
<font size="2">
<p>
	Datasets collected in association with the paper are linked below in the format of "<b>file:[file type]</b>".<br>
	You can read about file format requirements for Ymap by selecting the "Help" tab above.<br><br>

	Some datasets were retrieved from external databases and are linked below in the format of "<b>url:[database name, accession]</b>".<br>
	To download datasets from the NCBI-SRA database in FASTQ format, enter the experiment number for the sample into
	<a href="http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=search_seq_name">[this page]</a>, then select the accession number of interest.<br><br>
</p>
<p>
	<b>Figure 2 : WGseq, chr-end bias.</b><br>
	<ul>
		<li>YQ2, url:<a                     href="http://www.ebi.ac.uk/biosamples/sample/SAMEA1879786">[EMBL-EBI BiosSamples, SAMEA1879786]</a>
		    or url:<a                       href="http://www.ncbi.nlm.nih.gov/sra/ERX238051"          >[NCBI-SRA, ERX238051]</a></li>
	</ul>
	<b>Figure 3 : WGseq, GC-content bias.</b><br>
    <ul>
		<li>FH6, file:<a                    href="example_datasets/Fig3.FH6_forward.fastq.zip"        >[ZIP compressed FASTQ]</a>
		    and file:<a                     href="example_datasets/Fig3.FH6_reverse.fastq.zip"        >[ZIP compressed FASTQ]</a></li>
	</ul>
	<b>Figure 4 : ddRADseq.</b><br>
	<ul>
		<li>SC5314, file:<a                 href="example_datasets/Fig4.SC5314.txt"                   >[TXT]</a></li>
		<li>YJB11461, file:<a               href="example_datasets/Fig4.YJB11461.txt"                 >[TXT]</a></li>
		<li>YJB9442, file:<a                href="example_datasets/Fig4.YJB9442.txt"                  >[TXT]</a></li>
	</ul>
	<b>Figure 5 : WGseq</b><br>
	<ul>
		<li>SC5314, url:<a                  href="http://www.ncbi.nlm.nih.gov/sra/?term=SRR868699"    >[NCBI-SRA, SRR868699]</a></li>
		<li>FH1, file:<a                    href="example_datasets/Fig5.FH1_forward.fastq.zip"        >[ZIP compressed FASTQ]</a>
			and file:<a                     href="example_datasets/Fig5.FH1_reverse.fastq.zip"        >[ZIP compressed FASTQ]</a></li>
		<li>FH5, file:<a                    href="example_datasets/Fig5.FH5_forward.fastq.zip"        >[ZIP compressed FASTQ]</a>
		    and file:<a                     href="example_datasets/Fig5.FH5_reverse.fastq.zip"        >[ZIP compressed FASTQ]</a></li>
		<li>YJB12746, file:<a               href="example_datasets/Fig5.YJB12746_forward.fastq.zip"   >[ZIP compressed FASTQ]</a>
		    and file:<a                     href="example_datasets/Fig5.YJB12746_reverse.fastq.zip"   >[ZIP compressed FASTQ]</a></li>
	</ul>
	<b>Figure 6 : ddRADseq</b><br>
	<ul>
		<li>SC5314, file:<a                 href="example_datasets/Fig6.SC5314.txt"                   >[TXT]</a></li>
		<li>YJB12712 derivative #1, file:<a href="example_datasets/Fig6.YJB12712_derivative_1.txt"    >[TXT]</a></li>
		<li>YJB12712 derivative #2, file:<a href="example_datasets/Fig6.YJB12712_derivative_2.txt"    >[TXT]</a></li>
		<li>YJB12712 derivative #9, file:<a href="example_datasets/Fig6.YJB12712_derivative_9.txt"    >[TXT]</a></li>
	</ul>
	<b>Figure 8 : comparisons between array and sequence datasets</b><br>
	<ul>
		<li>YJB10490 array, file:<a         href="example_datasets/Fig8.YJB10490_array.xls"           >[XLS]</a></li>
		<li>YJB10490 WGseq, file:<a         href="example_datasets/Fig8.YJB10490_forward.fastq.zip"   >[ZIP compressed FASTQ]</a>
		    and file:<a                     href="example_datasets/Fig8.YJB10490_reverse.fastq.zip"   >[ZIP compressed FASTQ]</a></li>
		<li>YJB12229 array, file:<a         href="example_datasets/Fig8.YJB12229_array.xls"           >[XLS]</a></li>
		<li>YJB12229 ddRADseq, file:<a      href="example_datasets/Fig8.YJB12229.txt"                 >[TXT]</a></li>
		<li>Ss2 array, file:<a              href="example_datasets/Fig8.Ss2_array.xls"                >[XLS]</a></li>
		<li>YJB12353 WGseq, file:<a         href="example_datasets/Fig8.YJB12353_forward.fastq.zip"   >[ZIP compressed FASTQ</a>
		    and file:<a                     href="example_datasets/Fig8.YJB12353_reverse.fastq.zip"   >[ZIP compressed FASTQ]</a></li>
	</ul>
	<b>Figure 9 : LOH in clinical isolates</b><br>
	<ul>
		<li>SC5314, url:<a                  href="http://www.ebi.ac.uk/ena/data/view/SAMN02141741"    >[EMBL-EBI BioSamples, SAMN02141741]</a></li>
		<li>SC5314, file:<a                 href="example_datasets/Fig9.SC5314_forward.fastq.zip"     >[ZIP compressed FASTQ]</a>
		    and file:<a                     href="example_datasets/Fig9.SC5314_reverse.fastq.zip"     >[ZIP compressed FASTQ]</a></li>
		<li>SC5314, url:<a                  href="http://www.ncbi.nlm.nih.gov/sra/?term=SRR868699"    >[NCBI-SRA, SRR868699]</a></li>
		<li>FH1, file:<a                    href="example_datasets/Fig9.FH1_forward.fastq.zip"        >[ZIP compressed FASTQ]</a>
		    and file:<a                     href="example_datasets/Fig9.FH1_reverse.fastq.zip"        >[ZIP compressed FASTQ]</a></li>
		<li>ATCC200955, url:<a              href="http://www.ncbi.nlm.nih.gov/sra/?term=SAMN02140345" >[NCBI-SRA, SAMN02140345]</a></li>
		<li>ATCC10231, url:<a               href="http://www.ncbi.nlm.nih.gov/sra/?term=SAMN02140347" >[NCBI-SRA, SAMN02140347]</a></li>
		<li>YL1, url:<a                     href="http://www.ebi.ac.uk/biosamples/sample/SAMEA1879767">[EMBL-EBI BioSamples, SAMEA1879767]</a></li>
		<li>YQ2, url:<a                     href="http://www.ebi.ac.uk/biosamples/sample/SAMEA1879786">[EMBL-EBI BiosSamples, SAMEA1879786]</a>
		    or url:<a                       href="http://www.ncbi.nlm.nih.gov/sra/ERX238051"          >[NCBI-SRA, ERX238051]</a></li>
	</ul>
	<b>Figure 10 : FH series</b>
	<ul>
		<li>FH1, file:<a                    href="example_datasets/Fig10.FH1_forward.fastq.zip"       >[ZIP compressed FASTQ]</a>
		    and file:<a                     href="example_datasets/Fig10.FH1_reverse.fastq.zip"       >[ZIP compressed FASTQ]</a></li>
		<li>FH2, file:<a                    href="example_datasets/Fig10.FH2_forward.fastq.zip"       >[ZIP compressed FASTQ]</a>
		    and file:<a                     href="example_datasets/Fig10.FH2_reverse.fastq.zip"       >[ZIP compressed FASTQ]</a></li>
		<li>FH3, file:<a                    href="example_datasets/Fig10.FH3_forward.fastq.zip"       >[ZIP compressed FASTQ]</a>
		    and file:<a                     href="example_datasets/Fig10.FH3_reverse.fastq.zip"       >[ZIP compressed FASTQ]</a></li>
		<li>FH4, file:<a                    href="example_datasets/Fig10.FH4_forward.fastq.zip"       >[ZIP compressed FASTQ]</a>
		    and file:<a                     href="example_datasets/Fig10.FH4_reverse.fastq.zip"       >[ZIP compressed FASTQ]</a></li>
		<li>FH5, file:<a                    href="example_datasets/Fig10.FH5_forward.fastq.zip"       >[ZIP compressed FASTQ]</a>
		    and file:<a                     href="example_datasets/Fig10.FH5_reverse.fastq.zip"       >[ZIP compressed FASTQ]</a></li>
		<li>FH6, file:<a                    href="example_datasets/Fig10.FH6_forward.fastq.zip"       >[ZIP compressed FASTQ]</a>
		    and file:<a                     href="example_datasets/Fig10.FH6_reverse.fastq.zip"       >[ZIP compressed FASTQ]</a></li>
		<li>FH7, file:<a                    href="example_datasets/Fig10.FH7_forward.fastq.zip"       >[ZIP compressed FASTQ]</a>
		    and file:<a                     href="example_datasets/Fig10.FH7_reverse.fastq.zip"       >[ZIP compressed FASTQ]</a></li>
		<li>FH8, file:<a                    href="example_datasets/Fig10.FH8_forward.fastq.zip"       >[ZIP compressed FASTQ]</a>
		    and file:<a                     href="example_datasets/Fig10.FH8_reverse.fastq.zip"       >[ZIP compressed FASTQ]</a></li>
		<li>FH9, file:<a                    href="example_datasets/Fig10.FH9_forward.fastq.zip"       >[ZIP compressed FASTQ]</a>
		    and file:<a                     href="example_datasets/Fig10.FH9_reverse.fastq.zip"       >[ZIP compressed FASTQ]</a></li>
	</ul>
</p>
</font>
