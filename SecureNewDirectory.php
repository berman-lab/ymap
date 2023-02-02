<?php
//
// require_once SecureNewDirectory.php
//
function secureNewDirectory($dir) {
	// Generate 'index.php' file into each directory, to redirect to main page
	// to prevent users from exploring directory tree.
	$myfile = fopen($dir."index.php", "w");
	$txt = "<?php\nsession_start();\nerror_reporting(E_ALL);\nini_set('display_errors', 1);\nsession_destroy();\nheader('Location: ../');\n?>";
	fwrite($myfile,$txt);
	fclose($myfile);
}
?>
