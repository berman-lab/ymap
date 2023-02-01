<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	require_once 'POST_validation.php';
	ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		session_destroy();
		header('Location: .');
	}
	$user = $_SESSION['user'];

	// This script is intended to take information from file uploaders and then initiate the pipeline scripts to start processing.
	// This has been added to remove any server-side location information from being passed through client-side scripting.

	// validate POST strings.
	$dataFormat = sanitize_POST("dataFormat");
	$fileName   = sanitizeFile_POST("fileName");
	$genome     = sanitize_POST("genome");
	$project    = sanitize_POST("project");
	$key        = sanitize_POST("key");

	// Confirm if requested genome exists.
	$genome_dir = "users/".$user."/genomes/".$genome;
	if (!is_dir($genome_dir)) {
		// Genome doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: .');
	}

	// Confirm if requested project exists.
	$project_dir = "users/".$user."/projects/".$project;
	if (!is_dir($project_dir)) {
		// Project doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: .');
	}

	// set session variables.
	$_SESSION['dataFormat'] = $dataFormat;
	$_SESSION['fileName']   = $fileName;
	$_SESSION['genome']     = $genome;
	$_SESSION['project']    = $project;
	$_SESSION['key']        = $key;

	if ($project != "") {
		// initiate project processing.
		$conclusion_script = "";
		switch ($dataFormat) {
			case "SnpCghArray":
				$conclusion_script = "scripts_SnpCghArray/project.SnpCgh.install.php";
				break;
			case "WGseq_single":
				$conclusion_script = "scripts_seqModules/scripts_WGseq/project.single_WGseq.install_1.php";
				break;
			case "WGseq_paired":
				$conclusion_script = "scripts_seqModules/scripts_WGseq/project.paired_WGseq.install_1.php";
				break;
			case "ddRADseq_single":
				$conclusion_script = "scripts_seqModules/scripts_ddRADseq/project.single_ddRADseq.install_1.php";
				break;
			case "ddRADseq_paired":
				$conclusion_script = "scripts_seqModules/scripts_ddRADseq/project.paired_ddRADseq.install_1.php";
				break;
			case "RNAseq_single":
				$conclusion_script = "scripts_seqModules/scripts_RNAseq/project.single_RNAseq.install_1.php";
				break;
			case "RNAseq_paired":
				$conclusion_script = "scripts_seqModules/scripts_RNAseq/project.paired_RNAseq.install_1.php";
				break;
			case "IonExpressSeq_single":
				$conclusion_script = "scripts_seqModules/scripts_IonExpressSeq/project.single_IonExpressSeq.install_1.php";
				break;
			case "IonExpressSeq_paired":
				$conclusion_script = "scripts_seqModules/scripts_IonExpressSeq/project.paired_IonExpressSeq.install_1.php";
				break;
		}
	} else if ($genome != "") {
		// initiate genome processing.
		// ideally loaded into iframe id="Hidden_InstallNewGenome_Frame" defined in index.php
		$conclusion_script = "scripts_genomes/genome.install_1.php";
	} else {
		// No genome or project, should never happen: Force logout.
		session_destroy();
		header('Location: .');
	}

	// troubleshooting output
	//print "[upload_processer.php]\n";       print "user:        ".$user."\n";              print "project:     ".$project."\n";  print "genome:      ".$genome."\n";
	//print "data format: ".$dataFormat."\n"; print "script:      ".$conclusion_script."\n"; print "filename:    ".$fileName."\n"; print "key:         ".$key."\n";

	// Move to user directory
	chdir("users/".$user);

	// Open processing script.
	header("Location: ".$conclusion_script);
?>
