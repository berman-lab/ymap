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
	<BODY onload="UpdateProjectList()">
		<div id="hapmapCreationInformation"><p>
			<form action="scripts_seqModules/scripts_hapmaps/hapmap.install_1.php" method="post">
				<table><tr bgcolor="#CCFFCC"><td>
					<label for="hapmap">Hapmap Name : </label><input type="text" name="hapmap" id="hapmap">
				</td><td>
					Unique name for this hapmap.
				</td></tr><tr bgcolor="#CCCCFF"><td>
					<label for="genome">Reference Diploid Genome : </label><select name="genome" id="genome" onchange="UpdateProjectList()">
						<?php
						$genomesDir1    = "users/default/genomes/";
						$genomesDir2    = "users/".$user."/genomes/";
						$genomeFolders1 = array_diff(glob($genomesDir1."*"), array('..', '.'));
						$genomeFolders2 = array_diff(glob($genomesDir2."*"), array('..', '.'));
						foreach($genomeFolders1 as $key=>$folder) {   $genomeFolders1[$key] = str_replace($genomesDir1,"",$folder);   }
						foreach($genomeFolders2 as $key=>$folder) {   $genomeFolders2[$key] = str_replace($genomesDir2,"",$folder);   }
						$genomeFolders  = array_merge($genomeFolders1,$genomeFolders2);
						sort($genomeFolders);   //sort alphabetical.
						foreach ($genomeFolders as $key=>$genome) {
							echo "\n\t\t\t\t\t<option value='".$genome."'>".$genome."</option>";
						}
						?>
					</select><br>
				</td><td valign="top">
					Reference genome used to construct hapmap.
				</td></tr><tr bgcolor="#CCCCFF"><td>
					<label for="referencePloidy">Construct a hapmap from </label><select name="referencePloidy" id="referencePloidy" onchange="UpdateForm();">
					<option value="2">one diploid</option>
					<option value="1">two haploid</option>
					</select>
					<label for="referencePloidy"> strains.</label>
				</td><td>
					One heterozygous diploid or two haploids (or homozygous diploids)  can be used to define the SNPs of a hapmap.
				</td></tr><tr bgcolor="#CCFFCC"><td>
					<div id="hiddenFormSection1" style="display:inline">
<?php
// figure out which hapmaps have been defined for this species, if any.
$projectsDir1       = "users/default/projects/";
$projectsDir2       = "users/".$user."/projects/";
$projectFolders1    = array_diff(glob($projectsDir1."*"), array('..', '.'));
$projectFolders2    = array_diff(glob($projectsDir2."*"), array('..', '.'));
$projectFolders_raw = array_merge($projectFolders1,$projectFolders2);
// Go through each $projectFolder and look at 'genome.txt' and 'dataType.txt'; build javascript array of prejectName:genome:datatype triplets.
?>
						Parental strain : <select id="parent" name="parent"><option>[choose]</option></select>
<script type="text/javascript">
var projectGenomeDatatype_entries = [['project','genome','dataType']<?php
foreach ($projectFolders_raw as $key=>$folder) {
	$handle1         = fopen($folder."/genome.txt", "r");
	$genome_string   = trim(fgets($handle1));
	fclose($handle1);
	$handle2         = fopen($folder."/dataType.txt", "r");
	$dataType_string = trim(fgets($handle2));
	$dataType_string = explode(":",$dataType_string);
	$dataType_string = $dataType_string[0];
	fclose($handle2);
	$projectName     = $folder;
	$projectName     = str_replace($projectsDir1,"",$projectName);
	$projectName     = str_replace($projectsDir2,"",$projectName);
	echo ",['{$projectName}','{$genome_string}',{$dataType_string}]";
}
?>];

UpdateProjectList=function() {
	var selectedGenome   = document.getElementById("genome").value;     // grab genome.
	var selectedDatatype = 1
	var select           = document.getElementById("parent");     // grab select list.
	select.innerHTML     = '';
	for (var i = 1; i < projectGenomeDatatype_entries.length; i++) {
		var item = projectGenomeDatatype_entries[i];
		if (selectedGenome == item[1] && selectedDatatype == item[2]) {
			var el         = document.createElement("option");
			el.textContent = item[0];
			el.value       = item[0];
			select.appendChild(el);
		}
	}

	var select           = document.getElementById("child");
	select.innerHTML     = '';
	for (var i = 1; i < projectGenomeDatatype_entries.length; i++) {
		var item = projectGenomeDatatype_entries[i];
		if (selectedGenome == item[1] && selectedDatatype == item[2]) {
			var el         = document.createElement("option");
			el.textContent = item[0];
			el.value       = item[0];
			select.appendChild(el);
		}
	}

	var select           = document.getElementById("parentHaploid1");
	select.innerHTML     = '';
	for (var i = 1; i < projectGenomeDatatype_entries.length; i++) {
		var item = projectGenomeDatatype_entries[i];
		if (selectedGenome == item[1] && selectedDatatype == item[2]) {
			var el         = document.createElement("option");
			el.textContent = item[0];
			el.value       = item[0];
			select.appendChild(el);
		}
	}

	var select           = document.getElementById("parentHaploid2");
	select.innerHTML     = '';
	for (var i = 1; i < projectGenomeDatatype_entries.length; i++) {
		var item = projectGenomeDatatype_entries[i];
		if (selectedGenome == item[1] && selectedDatatype == item[2]) {
			var el         = document.createElement("option");
			el.textContent = item[0];
			el.value       = item[0];
			select.appendChild(el);
		}
	}
}
</script>
					</div>
				</td><td valign="top">
					<div id="hiddenFormSection2" style="display:inline">
						Only whole genome sequence datasets are used in constructing hapmaps.
					</div>
				</td></tr><tr bgcolor="#CCCCFF"><td valign="top">
					<div id="hiddenFormSection3" style="display:inline">
						First strain : <select id="child" name="child"><option>[choose]</option></select>
					</div<
				</td><td valign="top">
					<div id="hiddenFormSection4" style="display:inline">
						The first dataset used to construct the hapmap. (Others can be added later...)<br>
						Each strain used to construct the hapmap should have large loss of heterozygosity regions.
					</div>
				</td></tr><tr bgcolor="#CCFFCC"><td valign="top">
					<div id="hiddenFormSection5" style="display:none">
						Parental strain 1 : <select id="parentHaploid1" name="parentHaploid1"><option>[choose]</option></select>
					</div>
				</td><td valign="top">
					<div id="hiddenFormSection6" style="display:none">
						Only whole genome sequence datasets are used in constructing hapmaps.<br>
						This strain will form haplotype 'a'.
					</div>
				</td></tr><tr bgcolor="#CCCCFF"><td valign="top">
					<div id="hiddenFormSection7" style="display:none">
						Parental strain 2 : <select id="parentHaploid2" name="parentHaploid2"><option>[choose]</option></select>
					</div>
				</td><td valign="top">
					<div id="hiddenFormSection8" style="display:none">
						This strain will form haplotype 'b'.
					</div<
				</td></tr></table><br>
				<input type="submit" value="Create New hapmap">
			</form>
		</p></div>
	</body>
</html>

<script type="text/javascript">
UpdateForm=function() {
	if (document.getElementById("referencePloidy").value == 2) {
		document.getElementById("hiddenFormSection1").style.display  = 'inline';
		document.getElementById("hiddenFormSection2").style.display  = 'inline';
		document.getElementById("hiddenFormSection3").style.display  = 'inline';
		document.getElementById("hiddenFormSection4").style.display  = 'inline';
		document.getElementById("hiddenFormSection5").style.display  = 'none';
		document.getElementById("hiddenFormSection6").style.display  = 'none';
		document.getElementById("hiddenFormSection7").style.display  = 'none';
		document.getElementById("hiddenFormSection8").style.display  = 'none';
	} else { // haploid.
		document.getElementById("hiddenFormSection1").style.display  = 'none';
		document.getElementById("hiddenFormSection2").style.display  = 'none';
		document.getElementById("hiddenFormSection3").style.display  = 'none';
		document.getElementById("hiddenFormSection4").style.display  = 'none';
		document.getElementById("hiddenFormSection5").style.display  = 'inline';
		document.getElementById("hiddenFormSection6").style.display  = 'inline';
		document.getElementById("hiddenFormSection7").style.display  = 'inline';
		document.getElementById("hiddenFormSection8").style.display  = 'inline';
	}
}
</script>
