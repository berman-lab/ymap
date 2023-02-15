<?php
	session_start();
	if (isset($_SESSION['logged_on'])) { $user = $_SESSION['user']; } else { $user = 'default'; }
?>
<style type="text/css">
	html * {
		font-family: arial !important;
	}
	.tab {
		margin-left:    1cm;
	}
</style>
<font size='3'><b>User Manual</b></font><br>
<font size="2">
<br>
<a href="YMAP_User_Manual.docx">A detailed description of the user interface and tools available in YMAP.</a><br>
</font>
<br><br>
<font size='3'><b>Descriptions of datatypes and input file format requirements.</b></font><br>
<font size="2">
<p>
	The Yeast Mapping Analysis Pipeline (YMAP) is intended to simplify the processing of whole genome array and sequence
	datasets from yeast (or other small genome organisms) into simple to interpret figures illustrating large-scale genomic
	alterations.
</p>
<p>
	The current version (v1.0) of the pipeline is designed to process the following types of data:
	<ol>
		<li>SnpCgh Microarray.</li>
		<li>Whole genome next-generation-seq (WGseq).</li>
		<li>Double digest restriction site associated sequencing (ddRADseq).</li>
	</ol>
</p>
<hr width="50%">
<p>
	The input file format options for each data type are as follows :
</p>
<p>
	<b>1. SnpCgh microarray</b>
	<ul>
		<li>Tab-delimited text (*.tdt) with the following format (standard output from BlueFuse for Microarrays v3.6):
			<div class="tab"><font size='2'>
				46 header rows, followed by data rows.<br>
				1st column: probe ID/name.<br>
				4th column: probe channel 1 data.<br>
				5th column: probe channel 2 data.<br>
				6th column: probe channel ratio data.<br>
				7th column: probe channel log<sub>2</sub>ratio data.</font>
			</div>
		</li>
		<li>Unused columns in each data lines should contain some text.</li>
	</ul>
</p>
<p>
	<b>2. WGseq</b> and <b>3. ddRADseq</b>
	<ul>
		<li>Single-end reads as one file in raw FASTQ (*.fastq) or compressed (*.zip; *.fastq.gz) format.</li>
		<li>Paired-end reads as two files in raw FASTQ (*.fastq) or two compressed (*.zip; *.fastq.gz) format.</li>
		<li>Sequence Alignment/Map file in raw (*.sam) or compressed (*.bam) format.</li>
		<li>For examination of custom-filtered data, tab-delimited text (*.tdt) with the following format can be used :
			<div class="tab"><font size='2'>
				1st column: chromosome name (must match selected reference).<br>
				2nd column: chromosome bp coordinate.<br>
				3rd column: most common base-call at coordinate.<br>
				4th column: count of most common base-call at coordinate.<br>
				5th column: 2nd most common base-call at coordinate.<br>
				6th column: Count of 2nd most common base-call at coordinate.<br>
				7th column: 3rd most common base-call at coordinate.<br>
				8th column: Count of 3rd most common base-call at coordinate.<br>
				9th column: 4th most common base-call at coordinate.<br>
				10th column: Count of 4th most common base-call at coordinate.</font>
			</div>
		</li>
	</ul>
</p>
<p>
	<b>3. New Genomes</b>
	<ul>
		<li>New reference genomes in FASTA (*.fasta) format can be installed to use YMAP with data from other species.</li>
		<li>Genomes much smaller or larger than the 16 Mbase of Candida albicans may not work well. Please contact admin to discuss testing.</li>
	</ul>
</p>
<p>
	Uploaded files named with other file types will be discarded.
	<ul>
		<li>*.tdt</li>
		<li>*.fasta</li>
		<li>*.fastq</li>
		<li>*.fastq.zip</li>
		<li>*.fastq.gz</li>
		<li>*.sam</li>
		<li>*.bam</li>
	</ul>
</p>
<hr width="50%">
<p>
	Projects are listed with those needing data file upload first (highlighted in red), followed by those which are in-process (highlighted in yellow), and then those
	that are complete (highlighted in green). As projects transition from one category to the next, the color of their line will be updated. After a page refresh,
	position in the list will be updated.
</p>
<p>
	Haplotype maps can be constructed starting from one heterozygous parent or two homozygous parent datasets. When starting from a single heterozygous parent, additional
	child datasets with large-scale homozygoses are needed to construct the final map.
</p>
</font>
