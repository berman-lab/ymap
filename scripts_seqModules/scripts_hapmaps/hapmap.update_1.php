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
	require_once '../../POST_validation.php';
        ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		session_destroy();
		?><script type="text/javascript"> parent.location.reload(); </script><?php
	}

	// Load user string from session.
	$user       = $_SESSION['user'];

	// Validate input strings.
	$hapmap     = sanitize_POST("hapmap");
	$genome     = sanitize_POST("genome");
	$parent     = sanitize_POST("parent");
	$selectNext = sanitize_POST("selectNext");

	// Confirm if requested hapmap exists.
	$hapmap_dir = "../../users/".$user."/hapmaps/".$hapmap;
	if (!is_dir($hapmap_dir)) {
		// Hapmap doesn't exist, should never happen: Force logout.
		session_destroy();
		?><script type="text/javascript"> parent.location.reload(); </script><?php
	}

	// Confirm if requested genome exists.
	$genome_dir1 = "../../users/".$user."/genomes/".$genome;
	$genome_dir2 = "../../users/default/genomes/".$genome;
	if (!(is_dir($genome_dir1) || is_dir($genome_dir2))) {
		// Genome doesn't exist, should never happen: Force logout.
		session_destroy();
		?><script type="text/javascript"> parent.location.reload(); </script><?php
	}

	// Confirm if requested parent project exists.
	$parent_dir1 = "../../users/".$user."/projects/".$parent;
	$parent_dir2 = "../../users/default/projects/".$parent;
	if (!(is_dir($parent_dir1) || is_dir($parent_dir2))) {
		// Parent project doesn't exist, should never happen: Force logout.
		session_destroy();
		?><script type="text/javascript"> parent.location.reload(); </script><?php
	}

	// Confirm if requested next project exists.
	$selectNext_dir1 = "../../users/".$user."/projects/".$selectNext;
	$selectNext_dir2 = "../../users/default/projects/".$selectNext;
	if (!(is_dir($selectNext_dir1) || is_dir($selectNext_dir2))) {
		// next project doesn't exist, should never happen: Force logout.
		session_destroy();
		?><script type="text/javascript"> parent.location.reload(); </script><?php
	}

	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";

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
