<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
	"http://w3.org/TR/html4/loose.dtd">
<HEAD>
	<style type="text/css">
		body {font-family: arial;}
	</style>
	<meta charset="UTF-8">
</HEAD>
<BODY>
<?php
    session_start();
    if(!isset($_SESSION['logged_on'])){ ?> <script type="text/javascript"> parent.reload(); </script> <?php } else { $user = $_SESSION['user']; }
    require_once '../../constants.php';
    echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";
	error_reporting(E_ALL);
	ini_set('display_errors', 1);

	// define some characters which shouldn't be in strings.
	$bad_chars       = array(".", ",", "\\", "/", " ", "'", '"', "(", ")", "[", "]");

	// Process data transfered to this page.
	$user            = $_SESSION['user'];
	$hapmap          = str_replace($bad_chars,"",trim( filter_input(INPUT_POST, "hapmap",     FILTER_SANITIZE_STRING) ));
	$genome          = str_replace($bad_chars,"",trim( filter_input(INPUT_POST, "genome",     FILTER_SANITIZE_STRING) ));
	$parent          = str_replace($bad_chars,"",trim( filter_input(INPUT_POST, "parent",     FILTER_SANITIZE_STRING) ));
	$selectNext      = str_replace($bad_chars,"",trim( filter_input(INPUT_POST, "selectNext", FILTER_SANITIZE_STRING) ));

	// Make a form to generate a form to POST information to pass along to the next page in the process.
	echo "<script type=\"text/javascript\">\n";
	echo "\tvar autoSubmitForm = document.createElement('form');\n";
	echo "\t\tautoSubmitForm.setAttribute('method','post');\n";
	echo "\t\tautoSubmitForm.setAttribute('action','hapmap.update_2.php');\n";

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

	$project1 = $parent;
	$project2 = $selectNext;

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
	echo "\t\tdocument.body.appendChild(autoSubmitForm);";
	echo "\tautoSubmitForm.submit();\n";
	echo "</script>";
?>
</BODY>
</HTML>
