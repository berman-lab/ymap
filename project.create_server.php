<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		header('Location: user.login.php');
	}

	// Process POST data.
	$bad_chars              = array("~","@","#","$","%","^","&","*","(",")","+","=","|","{","}","<",">","?",".",",","\\","/","'",'"',"[","]","!");
	$projectName = trim(filter_input(INPUT_POST, "project", FILTER_SANITIZE_STRING));
	// removing unwanted characters
	$projectNameTrimmed     = str_replace($bad_chars,"",$projectName);
	// changing spaces to underlines
	$projectNameTrimmed     = str_replace(" ","_",$projectNameTrimmed);
	$ploidy                 = filter_input(INPUT_POST, "ploidy",                   FILTER_SANITIZE_STRING);
	$ploidyBase             = filter_input(INPUT_POST, "ploidyBase",               FILTER_SANITIZE_STRING);
	$dataFormat             = filter_input(INPUT_POST, "dataFormat",               FILTER_SANITIZE_STRING);
	$showAnnotations        = filter_input(INPUT_POST, "showAnnotations",          FILTER_SANITIZE_STRING);
	if ($dataFormat != "0") {		// we're dealing with ddRADseq, WGseq, or RNAseq.
		$readType            = filter_input(INPUT_POST, "readType",                 FILTER_SANITIZE_STRING);
		$performIndelRealign = filter_input(INPUT_POST, "indelrealign",             FILTER_SANITIZE_STRING);
		if ($performIndelRealign == "True") {
			$indelRealign = "1";
		} else {
			$indelRealign = "0";
		}
		$genome              = filter_input(INPUT_POST, "genome",                   FILTER_SANITIZE_STRING);
		$restrictionEnzymes  = filter_input(INPUT_POST, "selectRestrictionEnzymes", FILTER_SANITIZE_STRING);
		$hapmap              = filter_input(INPUT_POST, "selectHapmap",             FILTER_SANITIZE_STRING);
		$parent              = filter_input(INPUT_POST, "selectParent",             FILTER_SANITIZE_STRING);
		if (($parent == "none") || ($parent == "")) {
			$parent      = $projectNameTrimmed;  // no parent is used, so all calculations use the current project name as the parent.
		}
	}
	$manualLOH         = filter_input(INPUT_POST, "manualLOH",                    FILTER_SANITIZE_STRING);

	$user              = $_SESSION["user"];
	$dir1              = "users/".$user."/projects";
	$dir2              = "users/".$user."/projects/".$projectNameTrimmed;
	$dir3              = "users/default/projects/".$projectNameTrimmed;

	// Deals with accidental deletion of projects dir.
	if (!file_exists($dir1)){
		mkdir($dir1);
		chmod($dir1,0777);
	}

	if (file_exists($dir2) || file_exists($dir3)) {
		// Directory already exists
		echo "Project '".$projectName."' directory already exists.";
?>
	<html>
	<body>
	<script type="text/javascript">
	var el1 = parent.document.getElementById('Hidden_InstallNewDataset');
	el1.style.display = 'none';

	var el2 = parent.document.getElementById('panel_manageDataset_iframe').contentDocument.getElementById('name_error_comment');
	el2.style.visibility = 'visible';

	window.location = "project.create_window.php";
	</script>
	</body>
	</html>
<?php
	} else {
		// Create the project folder inside the user's projects directory
		mkdir($dir2);
		chmod($dir2,0777);

		// Generate 'name.txt' file containing:
		//      one line; name of genome.
		$outputName   = "users/".$user."/projects/".$projectNameTrimmed."/name.txt";
		$output       = fopen($outputName, 'w');
		fwrite($output, $projectName);
		fclose($output);

		$_SESSION['pending_install_project_count'] += 1;

		// Generate 'ploidy.txt' file.
		$fileName = "users/".$user."/projects/".$projectNameTrimmed."/ploidy.txt";
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
		$fileName = "users/".$user."/projects/".$projectNameTrimmed."/parent.txt";
		$file     = fopen($fileName, 'w');
		if ($dataFormat == "0") {
			fwrite($file, "none");
		} else {
			fwrite($file, $parent);
		}
		fclose($file);
		chmod($fileName,0644);

		// Generate 'dataFormat.txt' and 'dataBiases.txt' files.
		// dataFormat.txt file: #:#:# where 1st # indicates type of data, 2nd # indicates format of input data, & 3rd # indicates if indel-realignment should be done.
		// 1st #: 0=SnpCghArray; 1=WGseq; 2=ddRADseq; 3=RNAseq; 4=IonExpressSeq.
		// 2nd #: 0=single-end-reads FASTQ/ZIP/GZ; 1=paired-end-reads FASTQ/ZIP/GZ; 2=SAM/BAM; 3=TXT.
		// 3rd #: 0=False, no indel-realignment; 1=True, performe indel-realignment.
		$fileName1 = "users/".$user."/projects/".$projectNameTrimmed."/dataFormat.txt";
		$file1     = fopen($fileName1, 'w');
		$fileName2 = "users/".$user."/projects/".$projectNameTrimmed."/dataBiases.txt";
		$file2     = fopen($fileName2, 'w');
		if ($dataFormat == "0") { // SnpCghArray
			fwrite($file1, $dataFormat);
			$bias_GC     = filter_input(INPUT_POST, "0_bias2", FILTER_SANITIZE_STRING);
			$bias_end    = filter_input(INPUT_POST, "0_bias4", FILTER_SANITIZE_STRING);
			if (strcmp($bias_GC ,"") == 0) { $bias_GC  = "False"; }
			if (strcmp($bias_end,"") == 0) { $bias_end = "False"; }
			fwrite($file2,"False\n".$bias_GC."\nFalse\n".$bias_end);
		} else if ($dataFormat == "1") { // WGseq
			fwrite($file1, $dataFormat.":".$readType.":".$indelRealign);
			$bias_GC     = filter_input(INPUT_POST, "1_bias2", FILTER_SANITIZE_STRING);
			$bias_end    = filter_input(INPUT_POST, "1_bias4", FILTER_SANITIZE_STRING);
			if (strcmp($bias_GC ,"") == 0) { $bias_GC  = "False"; }
			if (strcmp($bias_end,"") == 0) { $bias_end = "False"; } else {$bias_GC  = "True"; }
			fwrite($file2,"False\n".$bias_GC."\nFalse\n".$bias_end);
		} else if ($dataFormat == "4") { // IonExpressSeq
			fwrite($file1, $dataFormat.":".$readType.":".$indelRealign);
			$bias_GC     = filter_input(INPUT_POST, "4_bias2", FILTER_SANITIZE_STRING);
			$bias_end    = filter_input(INPUT_POST, "4_bias4", FILTER_SANITIZE_STRING);
			if (strcmp($bias_GC ,"") == 0) { $bias_GC  = "False"; }
			if (strcmp($bias_end,"") == 0) { $bias_end = "False"; }
			fwrite($file2,"False\n".$bias_GC."\nFalse\n".$bias_end);
		} else if ($dataFormat == "2") { // ddRADseq
			fwrite($file1, $dataFormat.":".$readType.":".$indelRealign);
			$bias_length = filter_input(INPUT_POST, "2_bias1", FILTER_SANITIZE_STRING);
			$bias_GC     = filter_input(INPUT_POST, "2_bias2", FILTER_SANITIZE_STRING);
			$bias_end    = filter_input(INPUT_POST, "2_bias4", FILTER_SANITIZE_STRING);
			if (strcmp($bias_length,"") == 0) { $bias_length = "False"; }
			if (strcmp($bias_GC    ,"") == 0) { $bias_GC     = "False"; }
			if (strcmp($bias_end   ,"") == 0) { $bias_end    = "False"; }
			fwrite($file2,$bias_length."\n".$bias_GC."\nFalse\n".$bias_end);
		} else if ($dataFormat == "3") { // RNAseq
			fwrite($file1, $dataFormat.":".$readType.":".$indelRealign);
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

		// Generate 'restrictionEnzymes.txt' file, only for ddRADseq projects.
		if ($dataFormat == "2") { // ddRADseq
			$fileName = "users/".$user."/projects/".$projectNameTrimmed."/restrictionEnzymes.txt";
			$file     = fopen($fileName, 'w');
			fwrite($file, $restrictionEnzymes);
			fclose($file);
			chmod($fileName,0644);
		}

		// Generate 'snowAnnotations.txt' file.
		$fileName = "users/".$user."/projects/".$projectNameTrimmed."/showAnnotations.txt";
		$file     = fopen($fileName, 'w');
		fwrite($file, $showAnnotations);
		fclose($file);
		chmod($fileName,0644);

		// Generate 'genome.txt' file : containing genome used.
		//	1st line : (String) genome name.
		//	2nd line : (String) hapmap name.
		$fileName = "users/".$user."/projects/".$projectNameTrimmed."/genome.txt";
		$file     = fopen($fileName, 'w');
		if ($hapmap == "none") {
			fwrite($file, $genome);
		} else {
			fwrite($file, $genome."\n".$hapmap);
		}
		fclose($file);
		chmod($fileName,0644);

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
			$fileName = "users/".$user."/projects/".$projectNameTrimmed."/manualLOH.txt";
			$file     = fopen($fileName, 'w');
			fwrite($file, $manualLOH);
			fclose($file);
			chmod($fileName,0644);
		}
?>
	<html>
	<body>
	<script type="text/javascript">
	var el1 = parent.document.getElementById('panel_manageDataset_iframe').contentDocument.getElementById('newly_installed_list');
	el1.innerHTML += "<?php echo $_SESSION['pending_install_project_count']; ?>. <?php echo $projectName; ?><br>";

	var el2 = parent.document.getElementById('panel_manageDataset_iframe').contentDocument.getElementById('pending_comment');
	el2.style.visibility = 'visible';

	var el3 = parent.document.getElementById('panel_manageDataset_iframe').contentDocument.getElementById('name_error_comment');
        el3.style.visibility = 'hidden';

	var el4 = parent.document.getElementById('Hidden_InstallNewDataset');
	el4.style.display = 'none';

	window.location = "project.create_window.php";
	</script>
	</body>
	</html>
<?php
	}
?>
