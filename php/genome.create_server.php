<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.


	if(!isset($_SESSION['logged_on'])){
		header('Location: '.$GLOBALS['url'].'login.php');
	}

	$bad_chars      = array("~","@","#","$","%","^","&","*","(",")","+","=","|","{","}","<",">","?",".",",","\\","/"," ","'",'"',"[","]","!");
	$genomeName     = str_replace($bad_chars,"",trim( filter_input(INPUT_POST, "newGenomeName", FILTER_SANITIZE_STRING) ));
	$user           = $_SESSION['user'];
	$dir1           = $GLOBALS['directory']."users/".$user."/genomes";
	$dir2           = $GLOBALS['directory']."users/".$user."/genomes/".$genomeName;
	$dir3           = $GLOBALS['directory']."users/default/genomes/".$genomeName;

	// Deals with accidental deletion of genomes dir.
	if (!file_exists($dir1)){
		mkdir($dir1);
		chmod($dir1,0777);
	}

	if (file_exists($dir2) || file_exists($dir3)) {
		// Directory already exists
		echo "Genome '".$genomeName."' directory already exists.";
	} else {
		// Create the genome folder inside the user's genomes directory
		mkdir($dir2);
		chmod($dir2,0777);
		header('Location: '.$GLOBALS['url']);

		// Generate 'name.txt' file containing:
	    //      one line; name of genome.
	    fwrite($logOutput, "\tGenerating 'name.txt' file.\n");
	    $outputName       = $directory."users/".$user."/genomes/".$genomeName."/name.txt";
		$output       = fopen($outputName, 'w');
		fwrite($output, str_replace("_"," ",$genomeName));
	    fclose($output);
	}
?>
