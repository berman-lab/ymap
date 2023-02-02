<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	require_once 'POST_validation.php';
	require_once 'SecureNewDirectory.php';
	ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		session_destroy();
		header('Location: .');
	}

	// Load user string from session.
	$user            = $_SESSION['user'];

	// Validate input strings.
	$project         = sanitize_POST("project");
	$ploidy          = sanitizeFloat_POST("ploidy");
	$ploidyBase      = sanitizeFloat_POST("ploidyBase");
	$dataFormat      = sanitizeIntChar_POST("dataFormat");
	$showAnnotations = sanitizeIntChar_POST("showAnnotations");
	$manualLOH       = sanitizeTabbed_POST("manualLOH");

	// Validate additional inputs found when dealing with ddRADseq, WGseq, or RNAseq data types.
	if ($dataFormat != "0") {
		// Validate input strings.
		$readType            = sanitizeIntChar_POST("readType");
		$performIndelRealign = sanitize_POST("indelrealign");
			if ($performIndelRealign == "True") {	$indelRealign = "1";
			} else {				$indelRealign = "0";
			}
		$genome              = sanitize_POST("genome");
			// Confirm if requested genome exists.
			$genome_dir1 = "users/".$user."/genomes/".$genome;
			$genome_dir2 = "users/default/genomes/".$genome;
			if (!(is_dir($genome_dir1) || is_dir($genome_dir2))) {
				// Genome doesn't exist, should never happen: Force logout.
				session_destroy();
				header('Location: .');
			}
		$hapmap              = sanitize_POST("selectHapmap");
			if (($hapmap == "none") || ($hapmap == "")) {
				// no hapmap is used.
			} else {
				// Confirm if requested hapmap exists.
				$hapmap_dir1 = "users/".$user."/hapmaps/".$hapmap;
				$hapmap_dir2 = "users/default/hapmaps/".$hapmap;
				if (!(is_dir($hapmap_dir1) || is_dir($hapmap_dir2))) {
					// Hapmap doesn't exist, should never happen: Force logout.
					session_destroy();
					header('Location: .');
				}
			}
		$restrictionEnzymes  = sanitize_POST("selectRestrictionEnzymes");
		$parent              = sanitize_POST("selectParent");
			if (($parent == "none") || ($parent == "")) {
				// no parent is used, so all calculations use the current project name as the parent.
				$parent      = $project;
			} else {
				// Confirm if requested parent project exists.
				$parent_dir1 = "users/".$user."/projects/".$parent;
				$parent_dir2 = "users/default/projects/".$parent;
				if (!(is_dir($parent_dir1) || is_dir($parent_dir2))) {
					// Parent project doesn't exist, should never happen: Force logout.
					session_destroy();
					header('Location: .');
				}
			}
	}

	// Define some directories for later use.
	$projects_dir  = "users/".$user."/projects";
	$project_dir1  = "users/".$user."/projects/".$project;
	$project_dir2  = "users/default/projects/".$project;

	// Deals with accidental deletion of user/projects dir.
	if (!file_exists($projects_dir)){
		mkdir($projects_dir);
		secureNewDirectory($projects_dir);
		chmod($projects_dir,0777);
	}

	if (file_exists($project_dir1) || file_exists($project_dir2)) {
		//============================================
		// Project directory already exists, so exit.
		//--------------------------------------------
		echo "Project '".$project."' directory already exists.";
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
		//========================================================
		// Project directory doesn't exist, go about creating it.
		//--------------------------------------------------------

		// Create the project folder inside the user's projects directory
		mkdir($project_dir1);
		secureNewDirectory($project_dir1);
		chmod($project_dir1,0777);

		// Generate 'name.txt' file containing:
		//      one line; name of genome.
		$outputName   = $project_dir1."/name.txt";
		$output       = fopen($outputName, 'w');
		fwrite($output, $project);
		fclose($output);

		$_SESSION['pending_install_project_count'] += 1;

		// Generate 'ploidy.txt' file.
		$fileName = $project_dir1."/ploidy.txt";
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
		$fileName = $project_dir1."/parent.txt";
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
		$fileName1 = $project_dir1."/dataFormat.txt";
		$file1     = fopen($fileName1, 'w');
		$fileName2 = $project_dir1."/dataBiases.txt";
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
			$fileName = $project_dir1."/restrictionEnzymes.txt";
			$file     = fopen($fileName, 'w');
			fwrite($file, $restrictionEnzymes);
			fclose($file);
			chmod($fileName,0644);
		}

		// Generate 'snowAnnotations.txt' file.
		$fileName = $project_dir1."/showAnnotations.txt";
		$file     = fopen($fileName, 'w');
		fwrite($file, $showAnnotations);
		fclose($file);
		chmod($fileName,0644);

		// Generate 'genome.txt' file : containing genome used.
		//	1st line : (String) genome name.
		//	2nd line : (String) hapmap name.
		$fileName = $project_dir1."/genome.txt";
		$file     = fopen($fileName, 'w');
		if ($hapmap == "none") {
			fwrite($file, $genome);
		} else {
			fwrite($file, $genome."\n".$hapmap);
		}
		fclose($file);
		chmod($fileName,0644);

		// Generate 'manualLOH.txt' file : contains manual LOH annotation information.
		// one entry per line...  if input was provided.
		// tab-delimited channels.
		//    1. chrID
		//    2. startbp
		//    3. endbp
		//    4. R
		//    5. G
		//    6. B
		if (strlen($manualLOH) > 0) {
			$fileName = $project_dir1."/manualLOH.txt";
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
	el1.innerHTML += "<?php echo $_SESSION['pending_install_project_count']; ?>. <?php echo $project; ?><br>";

	var el2 = parent.document.getElementById('panel_manageDataset_iframe').contentDocument.getElementById('pending_comment');
	el2.style.visibility = 'visible';

	var el3 = parent.document.getElementById('panel_manageDataset_iframe').contentDocument.getElementById('name_error_comment');
        el3.style.visibility = 'hidden';

	var el4 = parent.document.getElementById('Hidden_InstallNewDataset');
	el4.style.display = 'none';

	window.location = "project.create_window.php";

	// Refresh "projectsShown" string;
	parent.update_projectsShown_after_new_project();
	</script>
	</body>
	</html>
<?php
	}
?>
