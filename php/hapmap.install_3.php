<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		header('Location: '.$GLOBALS['url'].'login.php');
	}

	$user            = $_SESSION['user'];
	$bad_chars       = array(".", ",", "\\", "/", " ");
	$hapmap          = str_replace($bad_chars,"_",trim( filter_input(INPUT_POST, "hapmap", FILTER_SANITIZE_STRING) ));
	$genome          = filter_input(INPUT_POST, "genome",             FILTER_SANITIZE_STRING);
	$referencePloidy = filter_input(INPUT_POST, "referencePloidy",    FILTER_SANITIZE_STRING);
	$project1        = filter_input(INPUT_POST, "project1",           FILTER_SANITIZE_STRING);
	$project2        = filter_input(INPUT_POST, "project2",           FILTER_SANITIZE_STRING);
	$colorA          = filter_input(INPUT_POST, "homolog_a_color",    FILTER_SANITIZE_STRING);
	$colorB          = filter_input(INPUT_POST, "homolog_b_color",    FILTER_SANITIZE_STRING);
	if ($referencePloidy == 2) {
		$hapmap_description = filter_input(INPUT_POST, "hapmap_description", FILTER_SANITIZE_STRING);
	}

	$dir1      = $GLOBALS['directory']."users/".$user."/hapmaps";
	$dir2      = $GLOBALS['directory']."users/".$user."/hapmaps/".$hapmap;
	$dir3      = $GLOBALS['directory']."users/default/hapmaps/".$hapmap;

// Create the hapmap folder inside the user's hapmaps directory
	if (file_exists($dir2)) {
	} else {
		mkdir($dir2);
		chmod($dir2,0777);
	}

// Initialize 'process_log.txt' file.
	$logOutputName = $directory."users/".$user."/hapmaps/".$hapmap."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "Running 'php/hapmap.install_2.php'.\n");

// Create 'colors.txt' file to contain colors used in haplotype figures.
	$handleName = $directory."users/".$user."/hapmaps/".$hapmap."/colors.txt";
	if (file_exists($handleName)) {
	} else {
		$handle     = fopen($handleName, 'w');
		fwrite($handle, $colorA."\n".$colorB);
		fclose($handle);
	}

// Create 'genome.txt' file to contain genome used in haplotype.
	$handleName = $directory."users/".$user."/hapmaps/".$hapmap."/genome.txt";
	if (file_exists($handleName)) {
	} else {
		$handle     = fopen($handleName, 'w');
		fwrite($handle, $genome);
		fclose($handle);
	}

// Create 'parent.txt' file to contain parent genome used in haplotype.
	$handleName = $directory."users/".$user."/hapmaps/".$hapmap."/parent.txt";
	if (file_exists($handleName)) {
	} else {
		$handle     = fopen($handleName, 'w');
		if ($referencePloidy == 2) {
			fwrite($handle, $project1);
		} else {
			fwrite($handle, $project1."\n".$project2);
		}
		fclose($handle);
	}

// Initialize 'process_log.txt' file.
	$logOutputName = $directory."users/".$user."/hapmaps/".$hapmap."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "Running 'php/hapmap.install_3.php'.\n");

// Initialize 'haplotype.txt' file to hold haplotype entry descriptions.
	if ($referencePloidy == 2) {
		$haplotypeFileName = $directory."users/".$user."/hapmaps/".$hapmap."/haplotypeMap.txt";
		$haplotypeFile     = fopen($haplotypeFileName, 'a');
		fwrite($haplotypeFile, $project1."\n".$project2."\n".$hapmap_description."\n");
		fclose($haplotypeFile);
	} else {
		// for haplotype map derived from two haploid strains, there will not be more entries after construction.
		// Because of this, 'parent.txt' will contain sufficient information to define haplotype map.
	}


// Generate 'working.txt' file to let pipeline know that processing is underway.
	$handleName      = $directory."users/".$user."/hapmaps/".$hapmap."/working.txt";
	$handle          = fopen($handleName, 'w');
	$startTimeString = date("Y-m-d H:i:s");
	fwrite($handle, $startTimeString);
	fclose($handle);

	fwrite($logOutput, "'sh/hapmap.install_3.php' completed.\n");
	fwrite($logOutput, "...........................................................\n");
	fclose($logOutput);
// Pass control over to a shell script ('sh/hapmap.install_4.sh') to continue processing and link with matlab.
//	$system_call_string = "sh ../sh/hapmap.install_4.sh ".$user." ".$hapmap." ".$directory." > /dev/null &";
	$system_call_string = "sh ../sh/hapmap.install_4.sh ".$user." ".$referencePloidy." ".$project1." ".$project2." ".$hapmap." ".$directory." > /dev/null &";
	system($system_call_string);

// Return to main page.
	header('Location: '.$GLOBALS['url']);
?>
