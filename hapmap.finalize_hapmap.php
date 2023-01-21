<?php
	session_start();
	error_reporting(E_ALL);
	require_once 'constants.php';
	ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		header('Location: user.login.php');
	}

	$user = $_SESSION['user'];
	$bad_chars = array("~","@","#","$","%","^","&","*","(",")","+","=","|","{","}","<",">","?",".",",","\\","/","'",'"',"[","]","!");
	$hapmap    = str_replace($bad_chars,"_",trim( filter_input(INPUT_POST, "hapmap", FILTER_SANITIZE_STRING) ));
	$key       = str_replace($bad_chars,"",trim(filter_input(INPUT_POST, "key",    FILTER_SANITIZE_STRING)));

	$dir1      = "users/".$user."/hapmaps";
	$dir2      = "users/".$user."/hapmaps/".$hapmap;
	$dir3      = "users/default/hapmaps/".$hapmap;

	// figure out what user the hapmap is installed under.
	$folder = "users/".$user."/hapmaps/".$hapmap."/";

	// Re-initialize 'process_log.txt' file.
	$logOutputName = "users/".$user."/hapmaps/".$hapmap."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "Log file restarted for hapmap finalization.\n");

	// Generate 'complete.txt' file to let the pipeline know that the haplotype map has been finalized.
	$handleName = "users/".$user."/hapmaps/".$hapmap."/complete.txt";
	$handle     = fopen($handleName, 'w');
	fwrite($handle, "complete");
	fclose($handle);

	// Run finalization python scripts to process Ymap-internal hapmap to final output versions.
	//	python [Ymap_root]/scripts_seqModules/scripts_hapmaps/process_hapmap.output_cleaned.py [raw_hapmap_file] > [cleaned_hapmap_file]
	//	python [Ymap_root]/scripts_seqModules/scripts_hapmaps/process_hapmap.output_errors.py  [raw_hapmap_file] > [errors_hapmap_file]

	fwrite($logOutput, "Haplotype map finalized.\n");
	fclose($logOutput);
?>
<script type="text/javascript">
	var ff = parent.parent.document.getElementById('panel_hapmap_iframe');
	ff.src = ff.src;
</script>
