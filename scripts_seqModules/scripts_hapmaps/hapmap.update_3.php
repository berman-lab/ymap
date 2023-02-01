<?php
	session_start();
	error_reporting(E_ALL);
        require_once '../../constants.php';
	require_once '../../POST_validation.php';
        ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		session_destroy();
		header('Location: ../../');
	}

	// Load user string from session.
	$user       = $_SESSION['user'];

	// Validate input strings.
	$hapmap             = validate_POST("hapmap");
	$genome             = validate_POST("genome");
	$project1           = validate_POST("project1");
	$project2           = validate_POST("project2");
	$colorA             = validate_POST("homolog_a_color");
	$colorB             = validate_POST("homolog_b_color");
	$hapmap_description = validate_POST("hapmap_description");

	// Confirm if requested hapmap exists.
	$hapmap_dir = "../../users/".$user."/hapmaps/".$hapmap;
	if !(is_dir($hapmap_dir)) {
		// Hapmap doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: ../../');
	}

	// Confirm if requested genome exists.
	$genome_dir1 = "../../users/".$user."/genomes/".$genome;
	$genome_dir2 = "../../users/default/genomes/".$genome;
	if !(is_dir($genome_dir1) || is_dir($genome_dir2)) {
		// Genome doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: ../../');
	}

	// Confirm if requested project1 project exists.
	$project1_dir1 = "../../users/".$user."/projects/".$project1;
	$project1_dir2 = "../../users/default/projects/".$project1;
	if !(is_dir($project1_dir1) || is_dir($project1_dir2)) {
		// Parent project doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: ../../');
	}

	// Confirm if requested project2 project exists.
	$project2_dir1 = "../../users/".$user."/projects/".$project2;
	$project2_dir2 = "../../users/default/projects/".$project2;
	if !(is_dir($project2_dir1) || is_dir($project2_dir2)) {
		// Project2 project doesn't exist, should never happen: Force logout.
		session_destroy();
		header('Location: ../../');
	}

	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";

	// Directories to use later.
	$currentPath = getcwd();

	// Initialize 'process_log.txt' file.
	$logOutputName = $hapmap_dir."/process_log.txt";
	$logOutput     = fopen($logOutputName, 'a');
	fwrite($logOutput, "Running 'scripts_seqModules/scripts_hapmaps/hapmap.update_3.php'.\n");
	fwrite($logOutput, "\tuser            = ".$user."\n");
	fwrite($logOutput, "\thapmap          = ".$hapmap."\n");
	fwrite($logOutput, "\tgenome          = ".$genome."\n");
	fwrite($logOutput, "\tproject1        = ".$project1."\n");
	fwrite($logOutput, "\tproject2        = ".$project2."\n");
	fwrite($logOutput, "\tcolorA          = ".$colorA."\n");
	fwrite($logOutput, "\tcolorB          = ".$colorB."\n");

	// Open 'haplotype.txt' file holding haplotype entry descriptions.
	// Add new haplotype entry descriptions.
	$haplotypeFileName = $hapmap_dir."/haplotypeMap.txt";
	$haplotypeFile     = fopen($haplotypeFileName, 'a');
	fwrite($haplotypeFile, $project1."\n".$project2."\n".$hapmap_description."\n");
	fclose($haplotypeFile);

	// Generate 'working.txt' file to let pipeline know that processing is underway.
	$handleName      = $hapmap_dir."/working.txt";
	$handle          = fopen($handleName, 'w');
	$startTimeString = date("Y-m-d H:i:s");
	fwrite($handle, $startTimeString);
	fclose($handle);

	fwrite($logOutput, "'scripts_seqModules/scripts_hapmaps/hapmap.update_3.php' completed.\n");
	fwrite($logOutput, "...........................................................\n");
	fwrite($logOutput, $currentPath."\n");
	fclose($logOutput);

	// Pass control over to a shell script ('scripts_seqModules/scripts_hapmaps/hapmap.update_4.sh') to continue processing and link with matlab.
	$system_call_string = "sh hapmap.update_4.sh ".$user." ".$project1." ".$project2." ".$hapmap." > /dev/null &";

	system($system_call_string);
?>
<script type="text/javascript">
	// Close user interface window involved in updating hapmap.
	var el3 = parent.document.getElementById('Hidden_AddToHapmap');
	el3.style.display = 'none';

	// Force reload of hapmap tab in main interface?
	var ff = parent.document.getElementById('panel_hapmap_iframe');
	ff.src = ff.src;
</script>
