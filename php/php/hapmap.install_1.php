<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		header('Location: '.$GLOBALS['url'].'user.login.php');
	}

	$bad_chars       = array(".", ",", "\\", "/", " ", "'", '"', "(", ")", "[", "]");
	$hapmap          = str_replace($bad_chars,"",trim( filter_input(INPUT_POST, "hapmap", FILTER_SANITIZE_STRING) ));
	$user            = $_SESSION['user'];
	$dir1            = $GLOBALS['directory']."users/".$user."/hapmaps";
	$dir2            = $GLOBALS['directory']."users/".$user."/hapmaps/".$hapmap;
	$dir3            = $GLOBALS['directory']."users/default/hapmaps/".$hapmap;

	$genome          = filter_input(INPUT_POST, "genome",          FILTER_SANITIZE_STRING);
	$referencePloidy = filter_input(INPUT_POST, "referencePloidy", FILTER_SANITIZE_STRING);
	// If ($referencePloidy == 2) use the following two variables.
	$parent          = filter_input(INPUT_POST, "parent",          FILTER_SANITIZE_STRING);
	$child           = filter_input(INPUT_POST, "child",           FILTER_SANITIZE_STRING);
	// If ($referencePloidy == 1) use the following two variables.
	$parentHaploid1  = filter_input(INPUT_POST, "parentHaploid1",  FILTER_SANITIZE_STRING);
	$parentHaploid2  = filter_input(INPUT_POST, "parentHaploid2",  FILTER_SANITIZE_STRING);

	// Deals with accidental deletion of hapmaps dir.
	if (!file_exists($dir1)){
		mkdir($dir1);
		chmod($dir1,0777);
	}

	if (file_exists($dir2) || file_exists($dir3)) {
		// Directory already exists
		echo "Hapmap '".$hapmap."' directory already exists.";
		exit;
	} else {
		// Make a form to generate a form to POST information to pass along to the next page in the process.
		echo "<script type=\"text/javascript\">\n";
		echo "\tvar autoSubmitForm = document.createElement('form');\n";
		echo "\t\tautoSubmitForm.setAttribute('method','post');\n";
		echo "\t\tautoSubmitForm.setAttribute('action','hapmap.install_2.php');\n";

		echo "\tvar input2 = document.createElement('input');\n";
		echo "\t\tinput2.setAttribute('type','hidden');\n";
		echo "\t\tinput2.setAttribute('name','hapmap');\n";
		echo "\t\tinput2.setAttribute('value','{$hapmap}');\n";
		echo "\t\tautoSubmitForm.appendChild(input2);\n";

		echo "\tvar input2 = document.createElement('input');\n";
		echo "\t\tinput2.setAttribute('type','hidden');\n";
		echo "\t\tinput2.setAttribute('name','genome');\n";
		echo "\t\tinput2.setAttribute('value','{$genome}');\n";
		echo "\t\tautoSubmitForm.appendChild(input2);\n";

		echo "\tvar input3 = document.createElement('input');\n";
		echo "\t\tinput3.setAttribute('type','hidden');\n";
		echo "\t\tinput3.setAttribute('name','referencePloidy');\n";
		echo "\t\tinput3.setAttribute('value','{$referencePloidy}');\n";
		echo "\t\tautoSubmitForm.appendChild(input3);\n";

		if ($referencePloidy == 2) {
			$project1 = $parent;
			$project2 = $child;
		} else {
			$project1 = $parentHaploid1;
			$project2 = $parentHaploid2;
		}

		echo "\tvar input4 = document.createElement('input');\n";
		echo "\t\tinput4.setAttribute('type','hidden');\n";
		echo "\t\tinput4.setAttribute('name','project1');\n";
		echo "\t\tinput4.setAttribute('value','{$project1}');\n";
		echo "\t\tautoSubmitForm.appendChild(input4);\n";

		echo "\tvar input5 = document.createElement('input');\n";
		echo "\t\tinput5.setAttribute('type','hidden');\n";
		echo "\t\tinput5.setAttribute('name','project2');\n";
		echo "\t\tinput5.setAttribute('value','{$project2}');\n";
		echo "\t\tautoSubmitForm.appendChild(input5);\n";

		echo "\tautoSubmitForm.submit();\n";
		echo "</script>";
	}
?>
