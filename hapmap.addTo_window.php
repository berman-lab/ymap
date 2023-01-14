<?php
	session_start();
	if(!isset($_SESSION['logged_on'])){ ?> <script type="text/javascript"> parent.reload(); </script> <?php }
	else {
		$user      = $_SESSION['user'];
		$key       = preg_replace('/\D/', '', $_GET['k']);  //Strip all non-numerical characters from string.
	}
	require_once 'constants.php';
	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";

	//============================================
	// Determine hapmap involved from $key value.
	//............................................
	$hapmapsDir    = "users/".$user."/hapmaps/";
	$hapmapFolders = array_diff(glob($hapmapsDir."*"), array('..', '.'));
	// Sort directories by date, newest first.
	array_multisort(array_map('filemtime', $hapmapFolders), SORT_DESC, $hapmapFolders);
	// Trim path from each folder string; find hapmap that matches the key received above.
	foreach($hapmapFolders as $key2=>$folder) {
		$hapmapFolders[$key2] = str_replace($hapmapsDir,"",$folder);
		if ($key2 == $key) { $hapmap = $hapmapFolders[$key2]; }
	}

	// Re-initialize 'process_log.txt' file.
	$logOutputName = "users/".$user."/hapmaps/".$hapmap."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "\n================================================\n");
	fwrite($logOutput, "Log file restarted for hapmap addition\n");
	fwrite($logOutput, "Running 'hapmap.addTo_1.php'.\n");

	//=====================================
	// Load genome and parent from hapmap.
	//.....................................
	// Defining directory location for later use.
	$folder = "users/".$user."/hapmaps/".$hapmap."/";
	// Load genome from 'hapmap/genome.txt'.
	$handle1 = fopen($folder."genome.txt", "r");
	$genome  = trim(fgets($handle1));
	fclose($handle1);
	// Load parent from 'hapmap/parent.txt'.
	$handle2 = fopen($folder."parent.txt", "r");
	$parent  = trim(fgets($handle2));
	fclose($handle2);
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
	<BODY onload="UpdateNextList()">
		<div id="hapmapCreationInformation"><p>
			<form action="scripts_seqModules/scripts_hapmaps/hapmap.update_1.php" method="post">
				<table><tr bgcolor="#CCFFCC"><td>
					<label for="hapmap">Hapmap Name : </label><input type="text" name="hapmap" id="hapmap" value="<?php echo $hapmap; ?>" readonly style="background-color: #EEFFEE;">
				</td><td>
					Unique name for this hapmap.
				</td></tr><tr bgcolor="#CCCCFF"><td>
					<label for="genome">Reference genome : </label><input type="text" name="genome" id="genome" value="<?php echo $genome; ?>" size="60" readonly style="background-color: #EEEEFF;">
				</td><td valign="top">
					Reference genome used to construct hapmap.
				</td></tr><tr bgcolor="#CCFFCC"><td>
					<label for="parent">Parental strain : </label><input type="text" name="parent" id="parent" value="<?php echo $parent; ?>" readonly style="background-color: #EEFFEE;">
				</td><td valign="top">
					Parent genome sequence used in the selected haplotype map.
				</td></tr><tr bgcolor="#CCCCFF"><td valign="top">
<?php
					// figure out which hapmaps have been defined for this species, if any.
					$projectsDir1       = "users/default/projects/";
					$projectsDir2       = "users/".$user."/projects/";
					$projectFolders1    = array_diff(glob($projectsDir1."*"), array('..', '.'));
					$projectFolders2    = array_diff(glob($projectsDir2."*"), array('..', '.'));
					$projectFolders_raw = array_merge($projectFolders1,$projectFolders2);
					// Go through each $projectFolder and look at 'genome.txt' and 'dataFormat.txt'; build javascript array of prejectName:genome:dataformat triplets.
					?>
					Next strain : <select id="selectNext" name="selectNext"><option>[choose]</option></select>
					<script type="text/javascript">
					var nextGenomeDataFormat_entries = [['next','genome','dataFormat']<?php
					foreach ($projectFolders_raw as $key=>$folder) {
						$handle2         = fopen($folder."/dataFormat.txt", "r");
						$dataFormat_string = trim(fgets($handle2));
						$dataFormat_string = explode(":",$dataFormat_string);
						$dataFormat_string = $dataFormat_string[0];
						fclose($handle2);

						// Exclude projects from unusable data types.
						if ($dataFormat_string == '0') {
							// 0 : array data is excluded from options.
						} elseif ($dataFormat_string == '1') {
							// 1 : WGseq data is usable.
							$handle1         = fopen($folder."/genome.txt", "r");
							$genome_string   = trim(fgets($handle1));
							fclose($handle1);

							$nextName        = $folder;
							$nextName        = str_replace($projectsDir1,"",$nextName);
							$nextName        = str_replace($projectsDir2,"",$nextName);
							echo ",['{$nextName}','{$genome_string}',{$dataFormat_string}]";
						} else {
							// 2 : ddRADseq data is unusable.
							// 3 : IonExpressSeq data is unusable.
							// 4 : RNAseq data is unusable.
						}
					}
					?>];

					// This function is run at time of page load, to update the list selection for the next strain to be added to hapmap.
					UpdateNextList=function() {
						var selectedGenome     = "<?php echo $genome; ?>";
						var selectedDataFormat = 1
						var select             = document.getElementById("selectNext");     // grab select list.
						select.innerHTML       = '';
						for (var i = 1; i < nextGenomeDataFormat_entries.length; i++) {
							var item = nextGenomeDataFormat_entries[i];
							if (selectedGenome == item[1] && selectedDataFormat == item[2]) {
								var el         = document.createElement("option");
								el.textContent = item[0];
								el.value       = item[0];
								select.appendChild(el);
							}
						}
					}
					</script>
				</td><td valign="top">
					The dataset being used to extend the previously defined hapmap. (Others can be added later...)<br>
					Each strain used to construct the hapmap should have large loss of heterozygosity regions.
				</td></tr></table><br>
				<input type="submit" value="Add entry to hapmap">
			</form>
		</p></div>
	</body>
</html>
