<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

	// Proces POST data.
	$bad_chars              = array("~","@","#","$","%","^","&","*","(",")","+","=","|","{","}","<",">","?",".",",","\\","/"," ","'",'"',"[","]","!");
	$name                   = str_replace($bad_chars,"",trim( filter_input(INPUT_POST, "project", FILTER_SANITIZE_STRING) )) ;
	$ploidy                 = filter_input(INPUT_POST, "ploidy",                   FILTER_SANITIZE_STRING);
	$ploidyBase             = filter_input(INPUT_POST, "ploidyBase",               FILTER_SANITIZE_STRING);
	$dataType               = filter_input(INPUT_POST, "dataType",                 FILTER_SANITIZE_STRING);
	$showAnnotations        = filter_input(INPUT_POST, "showAnnotations",          FILTER_SANITIZE_STRING);
	if ($dataType != "0") {		// we're dealing with ddRADseq, WGseq, or RNAseq.
		$readType           = filter_input(INPUT_POST, "readType",                 FILTER_SANITIZE_STRING);
		$genome             = filter_input(INPUT_POST, "genome",                   FILTER_SANITIZE_STRING);
		$restrictionEnzymes = filter_input(INPUT_POST, "selectRestrictionEnzymes", FILTER_SANITIZE_STRING);
		$hapmap             = filter_input(INPUT_POST, "selectHapmap",             FILTER_SANITIZE_STRING);
		$parent             = filter_input(INPUT_POST, "selectParent",             FILTER_SANITIZE_STRING);
		if (($parent == "none") || ($parent == "")) {
			$parent         = $name;  // no parent is used, so all calculations use the current project name as the parent.
		}
	}
	$manualLOH         = filter_input(INPUT_POST, "manualLOH",                    FILTER_SANITIZE_STRING);
	$user              = $_SESSION["user"];
	$dir1              = $GLOBALS["directory"]."users/".$user."/projects";
	$dir2              = $GLOBALS["directory"]."users/".$user."/projects/".$name;
	$directory         = $GLOBALS["directory"];

	// Deals with accidental deletion of projects dir.
	if (!file_exists($dir1)) {
		mkdir($dir1);
		chmod($dir1,0777);
	}
	if (file_exists($dir2)) {
		echo "Project '".$name."' directory already exists.<br>";
	} else {
	mkdir($dir2);
	chmod($dir2,0777);
//	header('Location: '.$GLOBALS['url'].'project.create.test.php');  //send browser back to main page.

	// Generate 'ploidy.txt' file.
		$fileName = $directory."users/".$user."/projects/".$name."/ploidy.txt";
		$file     = fopen($fileName, 'w');
		if (is_numeric($ploidy)) {
			fwrite($file, $ploidy."\n");
			if (is_numeric($ploidyBase)) {
				fwrite($file, $ploidyBase);
			} else {
				fwrite($file, "2.0");
			}
		} else {
			fwrite($file, "2.0\n");
			if (is_numeric($ploidy)) {
				fwrite($file, $ploidyBase);
			} else {
				fwrite($file, "2.0");
			}
		}
		fclose($file);
		chmod($fileName,0644);

	// Generate 'parent.txt' file.
		$fileName = $directory."users/".$user."/projects/".$name."/parent.txt";
		$file     = fopen($fileName, 'w');
		if ($dataType == "0") {
			fwrite($file, "none");
		} else {
			fwrite($file, $parent);
		}
		fclose($file);
		chmod($fileName,0644);

	// Generate 'dataType.txt' and 'dataBiases.txt' files.
		$fileName1 = $directory."users/".$user."/projects/".$name."/dataType.txt";
		$file1     = fopen($fileName1, 'w');
		$fileName2 = $directory."users/".$user."/projects/".$name."/dataBiases.txt";
		$file2     = fopen($fileName2, 'w');
		if ($dataType == "0") { // SnpCghArray
			fwrite($file1, $dataType);
			$bias_GC     = filter_input(INPUT_POST, "0_bias2", FILTER_SANITIZE_STRING);
			if (strcmp($bias_GC,"") == 0) { $bias_GC = "False"; }
			fwrite($file2,"False\n".$bias_GC."\nFalse\nFalse");
		} else if ($dataType == "1") { // WGseq
			fwrite($file1, $dataType.":".$readType);
			$bias_GC     = filter_input(INPUT_POST, "1_bias2", FILTER_SANITIZE_STRING);
			$bias_end    = filter_input(INPUT_POST, "1_bias4", FILTER_SANITIZE_STRING);
			if (strcmp($bias_GC ,"") == 0) { $bias_GC  = "False"; }
			if (strcmp($bias_end,"") == 0) { $bias_end = "False"; }
			fwrite($file2,"False\n".$bias_GC."\nFalse\n".$bias_end);
		} else if ($dataType == "4") { // IonExpress-seq
			fwrite($file1, $dataType.":".$readType);
			$bias_GC     = filter_input(INPUT_POST, "4_bias2", FILTER_SANITIZE_STRING);
			$bias_end    = filter_input(INPUT_POST, "4_bias4", FILTER_SANITIZE_STRING);
			if (strcmp($bias_GC ,"") == 0) { $bias_GC  = "False"; }
			if (strcmp($bias_end,"") == 0) { $bias_end = "False"; }
			fwrite($file2,"False\n".$bias_GC."\nFalse\n".$bias_end);
		} else if ($dataType == "2") { // ddRADseq
			fwrite($file1, $dataType.":".$readType);
			$bias_length = filter_input(INPUT_POST, "2_bias1", FILTER_SANITIZE_STRING);
			$bias_GC     = filter_input(INPUT_POST, "2_bias2", FILTER_SANITIZE_STRING);
			$bias_end    = filter_input(INPUT_POST, "2_bias4", FILTER_SANITIZE_STRING);
			if (strcmp($bias_length,"") == 0) { $bias_length = "False"; }
			if (strcmp($bias_GC    ,"") == 0) { $bias_GC     = "False"; }
			if (strcmp($bias_end   ,"") == 0) { $bias_end    = "False"; }
			fwrite($file2,$bias_length."\n".$bias_GC."\nFalse\n".$bias_end);
		} else if ($dataType == "3") { // RNAseq
			fwrite($file1, $dataType.":".$readType);
			$bias_length = filter_input(INPUT_POST, "3_bias1", FILTER_SANITIZE_STRING);
			$bias_GC     = filter_input(INPUT_POST, "3_bias2", FILTER_SANITIZE_STRING);
			$bias_end    = filter_input(INPUT_POST, "3_bias4", FILTER_SANITIZE_STRING);
			if (strcmp($bias_length,"") == 0) { $bias_length = "False"; }
			if (strcmp($bias_GC    ,"") == 0) { $bias_GC     = "False"; }
			if (strcmp($bias_end   ,"") == 0) { $bias_end    = "False"; }
			fwrite($file2,$bias_length."\n".$bias_GC."\nFalse\n".$bias_end);
		}
		fclose($file1);
		fclose($file2);
		chmod($fileName1,0644);
		chmod($fileName1,0644);

	// Generate 'restrictionEnzymes.txt' file, only fir ddRADseq projects.
		if ($dataType == "2") { // ddRADseq
			$fileName = $directory."users/".$user."/projects/".$name."/restrictionEnzymes.txt";
			$file     = fopen($fileName, 'w');
			fwrite($file, $restrictionEnzymes);
			fclose($file);
			chmod($fileName,0644);
		}

	// Generate 'snowAnnotations.txt' file.
		$fileName = $directory."users/".$user."/projects/".$name."/showAnnotations.txt";
		$file     = fopen($fileName, 'w');
		fwrite($file, $showAnnotations);
		fclose($file);
		chmod($fileName,0644);

	// Generate 'genome.txt' file : containing genome used.
	//	1st line : (String) genome name.
	//	2nd line : (String) hapmap name.
		if ($dataType != "0") {
			$fileName = $directory."users/".$user."/projects/".$name."/genome.txt";
			$file     = fopen($fileName, 'w');
			if ($hapmap == "none") {
				fwrite($file, $genome);
			} else {
				fwrite($file, $genome."\n".$hapmap);
			}
			fclose($file);
			chmod($fileName,0644);
		}

	// Generate 'manualLOH.txt' file : contains manual LOH annotation information.
	// one entry per pline...  if input was provided.
	// tab-delimited channels.
	//    1. chrID
	//    2. startbp
	//    3. endbp
	//    4. R
	//    5. G
	//    6. B
		if (strlen($manualLOH) > 0) {
			$fileName = $directory."users/".$user."/projects/".$name."/manualLOH.txt";
			$file     = fopen($fileName, 'w');
			fwrite($file, $manualLOH);
			fclose($file);
			chmod($fileName,0644);
		}

	}

	$_SESSION['pending_install_count'] += 1;
?>
<html>
<body>
<script type="text/javascript">
el = window.parent.document.getElementById('newly_installed_list');
el.innerHTML += "<?php echo $_SESSION['pending_install_count']; ?>. <?php echo $name; ?><br>";

el2 = window.parent.document.getElementById('pending_comment');
el2.style.visibility='visible';

parent.update_projectsShown_after_new_project();

window.location = "<?php echo $GLOBALS['url']; ?>project.create.test.php";
</script>
testing
</body>
</html>


<?php

?>
