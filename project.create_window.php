<?php
	session_start();
	if(!isset($_SESSION['logged_on'])){ ?> <script type="text/javascript"> parent.reload(); </script> <?php } else { $user = $_SESSION['user']; }
	require_once 'constants.php';
	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";
?>
<html lang="en">
	<HEAD>
		<style type="text/css">
			body {font-family: arial;}
			.tab {margin-left:   1cm;}
		</style>
		<meta http-equiv="content-type" content="text/html; charset=utf-8">
		<title>[Needs Title]</title>
	</HEAD>
	<BODY onload="UpdateHapmapList(); UpdateParentList()">
		<div id="loginControls"><p>
		</p></div>
		<div id="projectCreationInformation"><p>
			<form action="project.create_server.test.php" method="post">
				<table><tr bgcolor="#CCFFCC"><td>
					<label for="project">Dataset Name : </label><input type="text" name="project" id="project">
				</td><td>
					Unique name for this dataset.
				</td></tr><tr bgcolor="#CCCCFF"><td>
					<label for="ploidy">Ploidy of experiment : </label><input type="text" name="ploidy"  id="ploidy" value="2.0"><br>
				</td><td>
					A ploidy estimate for the strain being analyzed.
				</td></tr><tr bgcolor="#CCFFCC"><td>
					<label for="ploidy">Baseline ploidy : </label><input type="text" name="ploidyBase"  id="ploidyBase" value="2.0"><br>
				</td><td>
					The copy number to use as a baseline in drawing copy number variations.
				</td></tr><tr bgcolor="#CCCCFF"><td>
					<label for="showAnnotations">Generate figure with annotations?</label><select name="showAnnotations" id="showAnnotations">
						<option value="1">Yes</option>
						<option value="0">No</option>
					</select>
				</td><td>
					Genome annotations, such as rDNA locus, can be drawn at bottom of figures.
				</td></tr><tr bgcolor="#CCFFCC"><td>
					<label for="dataType">Data type : </label><select name="dataType" id="dataType" onchange="UpdateForm(); UpdateHapmap(); UpdateParentList()">
						<option value="0">SnpCgh microarray        </option>
						<option value="1">Whole genome NGS         </option>
						<option value="2">ddRADseq                 </option>
						<option value="3">RNAseq (testing)         </option>
						<option value="4">IonExpress-seq (testing) </option>
					<!--	<option value="5">RADseq           </option> --!>
					</select>
				</td><td>
					The type of data to be processed.
				</td></tr><tr bgcolor="#CCCCFF"><td valign="top">
					<div id="hiddenFormSection1" style="display:none">
						<label for="readType">Read type : </label><select name="readType" id="readType">
							<option value="0">single-end reads; FASTQ/ZIP/GZ file.</option>
							<option value="1">paired-end reads; FASTQ/ZIP/GZ files.</option>
							<option value="2">SAM/BAM file.</option>
							<option value="3">TXT file.</option>
							</select><br>
					</div>
				</td><td>
					<div id="hiddenFormSection2" style="display:none">
						Single-end or paired-end reads in FASTQ format can be compressed into ZIP or GZ archives or in SAM/BAM alignment files.<br>
						Tab-delimted TXT column format is described in 'About' tab of main page.
					</div>
				</td></tr><tr bgcolor="#CCFFCC"><td>
					<?php
					$genomesDir1       = "users/default/genomes/";
					$genomesDir2       = "users/".$user."/genomes/";
					$genomeFolders1    = array_diff(glob($genomesDir1."*"), array('..', '.'));
					$genomeFolders2    = array_diff(glob($genomesDir2."*"), array('..', '.'));
					foreach($genomeFolders1 as $key=>$folder) {   $genomeFolders1[$key] = str_replace($genomesDir1,"",$folder);   }
					foreach($genomeFolders2 as $key=>$folder) {   $genomeFolders2[$key] = str_replace($genomesDir2,"",$folder);   }
					$genomeFolders  = array_merge($genomeFolders1,$genomeFolders2);
					sort($genomeFolders);   //sort alphabetical.
					?>
					<div id="hiddenFormSection3" style="display:none">
						<label for="genome">Reference genome : </label><select name="genome" id="genome" onchange="UpdateHapmap(); UpdateHapmapList(); UpdateParentList()">
							<?php
							foreach ($genomeFolders as $key=>$genome) {
								echo "\n\t\t\t\t\t<option value='".$genome."'>".$genome."</option>";
							}
							?>
						</select><br>
					</div>
				</td><td valign="top">
					<div id="hiddenFormSection4" style="display:none">
					</div>
				</td></tr><tr bgcolor="#CCCCFF"><td>
					<?php
					// figure out which hapmaps have been defined for this species, if any.
					$hapmapsDir1       = $directory."users/default/hapmaps/";
					$hapmapsDir2       = $directory."users/".$user."/hapmaps/";
					$hapmapFolders1    = array_diff(glob($hapmapsDir1."*"), array('..', '.'));
					$hapmapFolders2    = array_diff(glob($hapmapsDir2."*"), array('..', '.'));
					$hapmapFolders_raw = array_merge($hapmapFolders1,$hapmapFolders2);
					// Go through each $hapmapFolder and look at 'genome.txt'; build javascript array of hapmapName:genome pairs.
					?>
					<div id="hiddenFormSection10" style="display:none">
						Restriction enzymes :
						<select id="selectRestrictionEnzymes" name="selectRestrictionEnzymes" onchange="UpdateParent();">
						<option value="MfeI_MboI">MfeI & MboI</option>
						<option value="BamHI_BclI">BamHI & BclI (testing)</option>
						</select>
					</div>
				</td><td valign="top">
					<div id="hiddenFormSection11" style="display:none">
						Analysis of ddRADseq data is limited to restriction fragments bound by both restriction enzymes.<br>
						If your restriction enzyme pair is not listed, you can request it be added by commenting to the admin in the "System" tab.
					</div>
				</td></tr><tr bgcolor="#CCCCFF"><td>
					<?php
					// figure out which hapmaps have been defined for this species, if any.
					$hapmapsDir1       = $directory."users/default/hapmaps/";
					$hapmapsDir2       = $directory."users/".$user."/hapmaps/";
					$hapmapFolders1    = array_diff(glob($hapmapsDir1."*"), array('..', '.'));
					$hapmapFolders2    = array_diff(glob($hapmapsDir2."*"), array('..', '.'));
					$hapmapFolders_raw = array_merge($hapmapFolders1,$hapmapFolders2);
					// Go through each $hapmapFolder and look at 'genome.txt'; build javascript array of hapmapName:genome pairs.
					?>
					<div id="hiddenFormSection5" style="display:none">
						Haplotype map : <select id="selectHapmap" name="selectHapmap" onchange="UpdateParent();"><option>[choose]</option></select>
						<script type="text/javascript">
						var hapmapGenome_entries = [['hapmap','genome']<?php
						foreach ($hapmapFolders_raw as $key=>$folder) {
							$handle        = fopen($folder."/genome.txt", "r");
							$genome_string = trim(fgets($handle));
							fclose($handle);
							$hapmapName    = $folder;
							$hapmapName    = str_replace($hapmapsDir1,"",$hapmapName);
							$hapmapName    = str_replace($hapmapsDir2,"",$hapmapName);
							echo ",['{$hapmapName}','{$genome_string}']";
						}
						?>];
						</script>
					</div>
				</td><td valign="top">
					<div id="hiddenFormSection6" style="display:none">
						A haplotype map defines the phasing of heterozygous SNPs across the genome and must be matched to the background of the experiment for informative results.
						SNP information from the hapmap will be used for SNP/LOH analsyses.
					</div>
				</td></tr><tr bgcolor="#CCFFCC"><td>
					<?php
					// figure out which hapmaps have been defined for this species, if any.
					$projectsDir1       = $directory."users/default/projects/";
					$projectsDir2       = $directory."users/".$user."/projects/";
					$projectFolders1    = array_diff(glob($projectsDir1."*"), array('..', '.'));
					$projectFolders2    = array_diff(glob($projectsDir2."*"), array('..', '.'));
					$projectFolders_raw = array_merge($projectFolders1,$projectFolders2);
					// Go through each $projectFolder and look at 'genome.txt' and 'dataType.txt'; build javascript array of prejectName:genome:datatype triplets.
					?>
					<div id="hiddenFormSection7" style="display:none">
						Parental strain : <select id="selectParent" name="selectParent"><option>[choose]</option></select>
						<script type="text/javascript">
						var parentGenomeDatatype_entries = [['parent','genome','dataType']<?php
						foreach ($projectFolders_raw as $key=>$folder) {
							$handle1         = fopen($folder."/genome.txt", "r");
							$genome_string   = trim(fgets($handle1));
							fclose($handle1);
					 		$handle2         = fopen($folder."/dataType.txt", "r");
							$dataType_string = trim(fgets($handle2));
							$dataType_string = explode(":",$dataType_string);
							$dataType_string = $dataType_string[0];
							fclose($handle2);
							$parentName      = $folder;
							$parentName      = str_replace($projectsDir1,"",$parentName);
							$parentName      = str_replace($projectsDir2,"",$parentName);
							echo ",['{$parentName}','{$genome_string}',{$dataType_string}]";
						}
						?>];
						</script>
					</div>
				</td><td valign="top">
					<div id="hiddenFormSection8a" style="display:none">
						This strain will act as the SNP distribution control.
					</div>
					<div id="hiddenFormSection8b" style="display:none">
						This strain will act as the CNV normalization control.
					</div>
				</td></tr><tr bgcolor="#CCFFCC"><td>
					<div id="hiddenFormSection9a" style="display:inline">
						<!-- SnpCgh array --!>
						<input type="checkbox"      name="0_bias2" value="True" checked>GC-content bias
					</div>
					<div id="hiddenFormSection9b" style="display:none">
						<!-- WGseq --!>
						<input type="checkbox"      name="1_bias2" value="True" checked>GC-content bias<br>
						<input type="checkbox"      name="1_bias4" value="True"        >chromosome-end bias
					</div>
					<div id="hiddenFormSection9c" style="display:none">
						<!-- ddRADseq --!>
						<input type="checkbox"      name="2_bias1" value="True" checked>fragment-length bias<br>
						<input type="checkbox"      name="2_bias2" value="True" checked>GC-content bias<br>
						<input type="checkbox"      name="2_bias4" value="True"        >chromosome-end bias
					</div>
					<div id="hiddenFormSection9d" style="display:none">
						<!-- RNAseq --!>
						<input type="checkbox"      name="3_bias1" value="True" checked>ORF-length bias<br>
						<input type="checkbox"      name="3_bias2" value="True" checked>GC-content bias<br>
						<input type="checkbox"      name="3_bias4" value="True"        >chromosome-end bias
					</div>
					<div id="hiddenFormSection9e" style="display:none">
						<!-- IonExpress-seq --!>
						<input type="checkbox"      name="4_bias2" value="True" checked>GC-content bias<br>
						<input type="checkbox"      name="4_bias4" value="True"        >chromosome-end bias
					</div>
				</td><td>
				Corrections applied.
				</td></tr></table><br>
				<input type="submit" value="Create New Dataset">
			</form>

			<script type="text/javascript">
			UpdateParent = function() {
				// if 'selectHapmap' isn't "[None defined]" then hide parental strain row.
				var selectedHapmap = document.getElementById("selectHapmap").value;
				if (selectedHapmap == 'none') {
					document.getElementById("hiddenFormSection7" ).style.display  = 'inline';
					document.getElementById("hiddenFormSection8a").style.display  = 'none';
					document.getElementById("hiddenFormSection8b").style.display  = 'inline';
				} else {
					if (document.getElementById("dataType").value == 2) {    // ddRADseq.
						document.getElementById("hiddenFormSection7" ).style.display  = 'inline';
						document.getElementById("hiddenFormSection8a").style.display  = 'none';
						document.getElementById("hiddenFormSection8b").style.display  = 'inline';
					} else {
						document.getElementById("hiddenFormSection7" ).style.display  = 'none';
						document.getElementById("hiddenFormSection8a").style.display  = 'none';
						document.getElementById("hiddenFormSection8b").style.display  = 'none';
					}
				}
			}
			UpdateHapmapList=function() {
				var selectedGenome = document.getElementById("genome").value;   // grab genome name.
				var select         = document.getElementById("selectHapmap");   // grab select list.
				select.innerHTML   = '';
				var el             = document.createElement("option");
				el.textContent     = '[None defined]';
				el.value           = 'none';
				select.appendChild(el);
				for (var i = 1; i < hapmapGenome_entries.length; i++) {
					var item = hapmapGenome_entries[i];
					if (selectedGenome == item[1]) {
						var el         = document.createElement("option");
						el.textContent = item[0];
						el.value       = item[0];
						select.appendChild(el);
					}
				}
			}
			UpdateParentList=function() {
				var selectedGenome   = document.getElementById("genome").value;   // grab genome name.
				var selectedDatatype = document.getElementById("dataType").value; // grab dataset type.
				var select           = document.getElementById("selectParent");   // grab select list.
				select.innerHTML     = '';
				var el               = document.createElement("option");
				el.textContent       = '[This strain is parental type.]';
				el.value             = 'none';
				select.appendChild(el);
				for (var i = 1; i < parentGenomeDatatype_entries.length; i++) {
					var item = parentGenomeDatatype_entries[i];
					if (selectedGenome == item[1] && selectedDatatype == item[2]) {
						var el         = document.createElement("option");
						el.textContent = item[0];
						el.value       = item[0];
						select.appendChild(el);
					}
				}
			}
			UpdateForm=function() {
				if (document.getElementById("dataType").value == 0) {			// SnpCgh Microarray.
					document.getElementById("hiddenFormSection1").style.display  = 'none';
					document.getElementById("hiddenFormSection2").style.display  = 'none';
					document.getElementById("hiddenFormSection3").style.display  = 'none';
					document.getElementById("hiddenFormSection4").style.display  = 'none';
					document.getElementById("hiddenFormSection5").style.display  = 'none';
					document.getElementById("hiddenFormSection6").style.display  = 'none';
					document.getElementById("hiddenFormSection7").style.display  = 'none';
					document.getElementById("hiddenFormSection9a").style.display = 'inline';
					document.getElementById("hiddenFormSection9b").style.display = 'none';
					document.getElementById("hiddenFormSection9c").style.display = 'none';
					document.getElementById("hiddenFormSection9d").style.display = 'none';
					document.getElementById("hiddenFormSection9e").style.display = 'none';
					document.getElementById("hiddenFormSection10").style.display = 'none';
					document.getElementById("hiddenFormSection11").style.display = 'none';
				} else {														// WGseq or ddRADseq.
					document.getElementById("hiddenFormSection1").style.display  = 'inline';
					document.getElementById("hiddenFormSection2").style.display  = 'inline';
					document.getElementById("hiddenFormSection3").style.display  = 'inline';
					document.getElementById("hiddenFormSection4").style.display  = 'inline';
					document.getElementById("hiddenFormSection5").style.display  = 'inline';
					document.getElementById("hiddenFormSection6").style.display  = 'inline';
					document.getElementById("hiddenFormSection7").style.display  = 'inline';
					document.getElementById("hiddenFormSection10").style.display = 'none';
					document.getElementById("hiddenFormSection11").style.display = 'none';
					if (document.getElementById("dataType").value == 1) { // WGseq
						document.getElementById("hiddenFormSection9a").style.display = 'none';
						document.getElementById("hiddenFormSection9b").style.display = 'inline';
						document.getElementById("hiddenFormSection9c").style.display = 'none';
						document.getElementById("hiddenFormSection9d").style.display = 'none';
						document.getElementById("hiddenFormSection9e").style.display = 'none';
					} else if (document.getElementById("dataType").value == 2) { // ddRADseq
						document.getElementById("hiddenFormSection9a").style.display = 'none';
						document.getElementById("hiddenFormSection9b").style.display = 'none';
						document.getElementById("hiddenFormSection9c").style.display = 'inline';
						document.getElementById("hiddenFormSection9d").style.display = 'none';
						document.getElementById("hiddenFormSection9e").style.display = 'none';
						document.getElementById("hiddenFormSection10").style.display = 'inline';
						document.getElementById("hiddenFormSection11").style.display = 'inline';
					} else if (document.getElementById("dataType").value == 3) { // RNAseq (tsting)
						document.getElementById("hiddenFormSection9a").style.display = 'none';
						document.getElementById("hiddenFormSection9b").style.display = 'none';
						document.getElementById("hiddenFormSection9c").style.display = 'none';
						document.getElementById("hiddenFormSection9d").style.display = 'inline';
						document.getElementById("hiddenFormSection9e").style.display = 'none';
					} else { // IonExpress (testing)
						document.getElementById("hiddenFormSection9a").style.display = 'none';
						document.getElementById("hiddenFormSection9b").style.display = 'none';
						document.getElementById("hiddenFormSection9c").style.display = 'none';
						document.getElementById("hiddenFormSection9d").style.display = 'none';
						document.getElementById("hiddenFormSection9e").style.display = 'inline';
					}
				}
			}
			UpdateHapmap=function() {
				if (document.getElementById("dataType").value == 0) {			// SnpCgh microarray.
					document.getElementById("hiddenFormSection8a").style.display = 'none';
					document.getElementById("hiddenFormSection8b").style.display = 'none';
				} else if (document.getElementById("dataType").value == 2) {	// ddRADseq.
					document.getElementById("hiddenFormSection8a").style.display = 'none';
					document.getElementById("hiddenFormSection8b").style.display = 'inline';
				} else {														// WGseq
					document.getElementById("hiddenFormSection8a").style.display = 'inline';
					document.getElementById("hiddenFormSection8b").style.display = 'none';
				}
			}
			</script>
		</p></div>
	</body>
</html>
