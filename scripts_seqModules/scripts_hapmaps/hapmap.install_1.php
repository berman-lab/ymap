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
	error_reporting(E_ALL);
        require_once '../../constants.php';
        ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		session_destroy();
		?> <script type="text/javascript"> parent.reload(); </script> <?php
	}

	// Load user string from session.
	$user       = $_SESSION['user'];

	// Validate hapmap input string.
	$hapmap      = trim(filter_input(INPUT_POST, "hapmap", FILTER_SANITIZE_STRING));	// strip out any html tags.
	$hapmap      = str_replace(" ","_",$hapmap);					// convert any spaces to underlines.
	$hapmap      = preg_replace("/[\s\W]+/", "", $hapmap);				// remove everything but alphanumeric characters and underlines.
	$hapmap_dir1 = "../../users/".$user."/hapmaps/".$hapmap;
	$hapmap_dir2 = "../../users/default/hapmaps/".$hapmap;

	// Validate genome input string.
	$genome     = trim(filter_input(INPUT_POST, "genome", FILTER_SANITIZE_STRING));
	$genome     = str_replace(" ","_",$genome);
	$genome     = preg_replace("/[\s\W]+/", "", $genome);
	// Confirm if requested genome exists.
	$genome_dir1 = "../../users/".$user."/genomes/".$genome;
	$genome_dir2 = "../../users/default/genomes/".$genome;
	if !(is_dir($genome_dir1) || is_dir($genome_dir2)) {
		// Genome doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: user.login.php');
	}

	// Validate ploidy input string.
	$referencePloidy = trim(filter_input(INPUT_POST, "referencePloidy", FILTER_SANITIZE_STRING));	// strip out any html tags.
	$referencePloidy = preg_replace("/[^\d\.]+/", "", $ploidy);					// remove everything but numerals and period.

	// Validate parent input string.
	$parent      = trim(filter_input(INPUT_POST, "parent", FILTER_SANITIZE_STRING));
	$parent      = str_replace(" ","_",$parent);
	$parent      = preg_replace("/[\s\W]+/", "", $parent);
	// Confirm if requested project1 project exists.
	$parent_dir1 = "../../users/".$user."/projects/".$parent;
	$parent_dir2 = "../../users/default/projects/".$parent;
	if !(is_dir($parent_dir1) || is_dir($parent_dir2)) {
		// Parent project doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: user.login.php');
	}

	// Validate child input string.
	$child      = trim(filter_input(INPUT_POST, "child", FILTER_SANITIZE_STRING));
	$child      = str_replace(" ","_",$child);
	$child      = preg_replace("/[\s\W]+/", "", $child);
	// Confirm if requested project1 project exists.
	$child_dir1 = "../../users/".$user."/projects/".$child;
	$child_dir2 = "../../users/default/projects/".$child;
	if !(is_dir($child_dir1) || is_dir($child_dir2)) {
		// Child project doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: user.login.php');
	}

	// Validate parentHaploid1 input string.
	$parentHaploid1      = trim(filter_input(INPUT_POST, "parentHaploid1", FILTER_SANITIZE_STRING));
	$parentHaploid1      = str_replace(" ","_",$parentHaploid1);
	$parentHaploid1      = preg_replace("/[\s\W]+/", "", $parentHaploid1);
	// Confirm if requested project1 project exists.
	$parentHaploid1_dir1 = "../../users/".$user."/projects/".$parentHaploid1;
	$parentHaploid1_dir2 = "../../users/default/projects/".$parentHaploid1;
	if !(is_dir($parentHaploid1_dir1) || is_dir($parentHaploid1_dir2)) {
		// parentHaploid1 project doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: user.login.php');
	}

	// Validate parentHaploid2 input string.
	$parentHaploid2      = trim(filter_input(INPUT_POST, "parentHaploid2", FILTER_SANITIZE_STRING));
	$parentHaploid2      = str_replace(" ","_",$parentHaploid2);
	$parentHaploid2      = preg_replace("/[\s\W]+/", "", $parentHaploid2);
	// Confirm if requested project1 project exists.
	$parentHaploid2_dir1 = "../../users/".$user."/projects/".$parentHaploid2;
	$parentHaploid2_dir2 = "../../users/default/projects/".$parentHaploid2;
	if !(is_dir($parentHaploid2_dir1) || is_dir($parentHaploid2_dir2)) {
		// parentHaploid1 project doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: user.login.php');
	}

	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";

	// Deals with accidental deletion of user/hapmaps dir.
	$dir1            = "../../users/".$user."/hapmaps";
	if (!file_exists($dir1)){
		mkdir($dir1);
		chmod($dir1,0777);
	}

	if (file_exists($hapmap_dir1) || file_exists($hapmap_dir2)) {
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
		echo "\t\tdocument.body.appendChild(autoSubmitForm);";
		echo "\tautoSubmitForm.submit();\n";
		echo "</script>";
	}
?>
</BODY>
</HTML>
