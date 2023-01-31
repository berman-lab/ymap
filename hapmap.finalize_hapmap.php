<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		session_destroy();
		header('Location: user.login.php');
	}

	// Load user string from session.
	$user   = $_SESSION['user'];

	// Sanitize input strings.
	$hapmap = trim(filter_input(INPUT_POST, "hapmap", FILTER_SANITIZE_STRING));	// strip out any html tags.
	$hapmap = str_replace(" ","_",$hapmap);						// convert any spaces to underlines.
	$hapmap = preg_replace("/[\s\W]+/", "", $hapmap);				// remove everything but alphanumeric characters and underlines.
	$key    = trim(filter_input(INPUT_POST, "key", FILTER_SANITIZE_STRING));
	$key    = str_replace(" ","_",$key);
	$key    = preg_replace("/[\s\W]+/", "", $key);

	// Construct local dir strings.
	$hapmap_dir = "users/".$user."/hapmaps/".$hapmap;

	// Confirm if requested hapmap exists.
	if (!is_dir($hapmap_dir)) {
		// Hapmap doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: user.login.php');
        }

	// Re-initialize 'process_log.txt' file.
	$logOutputName = $hapmap_dir."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "Log file restarted for hapmap finalization.\n");

	// Generate 'complete.txt' file to let the pipeline know that the haplotype map has been finalized.
	$handleName = $hapmap_dir."/complete.txt";
	$handle     = fopen($handleName, 'w');
	fwrite($handle, "complete");
	fclose($handle);

	// log output.
	fwrite($logOutput, "Haplotype map finalized.\n");
	fclose($logOutput);
?>
<script type="text/javascript">
	var ff = parent.parent.document.getElementById('panel_hapmap_iframe');
	ff.src = ff.src;
</script>
