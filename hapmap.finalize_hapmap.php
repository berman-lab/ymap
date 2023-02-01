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

	// Load user string from session.
	$user   = $_SESSION['user'];

	// Sanitize input strings.
	$hapmap = sanitize_POST("hapmap");
	$key    = sanitize_POST("key");

	// Construct local dir strings.
	$hapmap_dir = "users/".$user."/hapmaps/".$hapmap;

	// Confirm if requested hapmap exists.
	if (!is_dir($hapmap_dir)) {
		// Hapmap doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: .');
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
