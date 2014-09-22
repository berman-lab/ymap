<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		header('Location: '.$GLOBALS['url'].'user.login.php');
	}

	$bad_chars = array(".", ",", "\\", "/", " ");
	$hapmap    = str_replace($bad_chars,"_",trim( filter_input(INPUT_POST, "hapmap", FILTER_SANITIZE_STRING) ));
	$user      = $_SESSION['user'];
	$genome    = filter_input(INPUT_POST, "genome",      FILTER_SANITIZE_STRING);
	$parent    = filter_input(INPUT_POST, "parent",      FILTER_SANITIZE_STRING);
	$child     = filter_input(INPUT_POST, "selectNext",  FILTER_SANITIZE_STRING);

	// Re-initialize 'process_log.txt' file.
	$logOutputName = $directory."users/".$user."/hapmaps/".$hapmap."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "Log file restarted for hapmap addition\n");
	fwrite($logOutput, "Running 'php/hapmap.addTo_1.php'.\n");

	// Make a form to generate a form to POST information to pass along to the next page in the process.
	echo "<script type=\"text/javascript\">\n";
	echo "\tvar autoSubmitForm = document.createElement('form');\n";
	echo "\t\tautoSubmitForm.setAttribute('method','post');\n";
	echo "\t\tautoSubmitForm.setAttribute('action','hapmap.install_2.php');\n";
	echo "\tvar input1 = document.createElement('input');\n";
	echo "\t\tinput1.setAttribute('type','hidden');\n";
	echo "\t\tinput1.setAttribute('name','hapmap');\n";
	echo "\t\tinput1.setAttribute('value','".$hapmap."');\n";
	echo "\t\tautoSubmitForm.appendChild(input1);\n";
	echo "\tvar input2 = document.createElement('input');\n";
	echo "\t\tinput2.setAttribute('type','hidden');\n";
	echo "\t\tinput2.setAttribute('name','genome');\n";
	echo "\t\tinput2.setAttribute('value','".$genome."');\n";
	echo "\t\tautoSubmitForm.appendChild(input2);\n";
	echo "\tvar input3 = document.createElement('input');\n";
	echo "\t\tinput3.setAttribute('type','hidden');\n";
	echo "\t\tinput3.setAttribute('name','project1');\n";
	echo "\t\tinput3.setAttribute('value','".$parent."');\n";
	echo "\t\tautoSubmitForm.appendChild(input3);\n";
	echo "\tvar input4 = document.createElement('input');\n";
	echo "\t\tinput4.setAttribute('type','hidden');\n";
	echo "\t\tinput4.setAttribute('name','project2');\n";
	echo "\t\tinput4.setAttribute('value','".$child."');\n";
	echo "\t\tautoSubmitForm.appendChild(input4);\n";
	echo "\tvar input5 = document.createElement('input');\n";
	echo "\t\tinput5.setAttribute('type','hidden');\n";
	echo "\t\tinput5.setAttribute('name','referencePloidy');\n";
	echo "\t\tinput5.setAttribute('value','2');\n";
	echo "\t\tautoSubmitForm.appendChild(input5);\n";
	echo "\tautoSubmitForm.submit();\n";
	echo "</script>";


	fwrite($logOutput, "'php/hapmap.addTo_1.php' completed.\n");
	fclose($logOutput);
?>
