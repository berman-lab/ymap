<?php
	session_start();
	error_reporting(E_ALL);
        require_once '../../constants.php';
	require_once '../../POST_validation.php';
	require_once '../../SecureNewDirectory.php';
        ini_set('display_errors', 1);

	// If the user is not logged on, redirect to login page.
	if(!isset($_SESSION['logged_on'])){
		session_destroy();
		?><script type="text/javascript"> parent.location.reload(); </script><?php
	}

	// Load user string from session.
	$user   = $_SESSION['user'];

	// Validate input strings.
	$hapmap          = sanitize_POST("hapmap");
	$genome          = sanitize_POST("genome");
	$project1        = sanitize_POST("project1");
	$project2        = sanitize_POST("project2");
	$colorA          = sanitizeColor_POST("homolog_a_color");
	$colorB          = sanitizeColor_POST("homolog_b_color");
	$referencePloidy = (float)sanitizeFloat_POST("referencePloidy");
	if ($referencePloidy == 2) {
		// Validate hapmap description string.
		$hapmap_description = sanitizeHapmap_POST("hapmap_description");
	}

	// Define hapmap directory.
	$hapmap_dir = "../../users/".$user."/hapmaps/".$hapmap;

	// Confirm if requested genome exists.
	$genome_dir1 = "../../users/".$user."/genomes/".$genome;
	$genome_dir2 = "../../users/default/genomes/".$genome;
	if (!(is_dir($genome_dir1) || is_dir($genome_dir2))) {
		// Genome doesn't exist, should never happen: Force logout.
		session_destroy();
		?><script type="text/javascript"> parent.location.reload(); </script><?php
	}

	// Confirm if requested projects exists.
	$project1_dir1   = "../../users/".$user."/projects/".$project1;
	$project1_dir2   = "../../users/default/projects/".$project1;
	$project2_dir1   = "../../users/".$user."/projects/".$project2;
	$project2_dir2   = "../../users/default/projects/".$project2;
	if (!(is_dir($project1_dir1) || is_dir($project1_dir2) || is_dir($project2_dir1) || is_dir($project2_dir2))) {
		// A project doesn't exist, should never happen: Force logout.
		session_destroy();
		?><script type="text/javascript"> parent.location.reload(); </script><?php
	}
	$currentPath = getcwd();

        if (file_exists($hapmap_dir)) {
		//============================================
		// Hapmap directory already exists, so exit.
		//--------------------------------------------
		echo "Hapmap '".$hapmap."' directory already exists.";
?>
        <html>
        <body>
        <script type="text/javascript">
        var el1 = parent.document.getElementById('Hidden_InstallNewDataset');
        el1.style.display = 'none';

        var el2 = parent.document.getElementById('panel_manageDataset_iframe').contentDocument.getElementById('name_error_comment');
        el2.style.visibility = 'visible';

        window.location = "project.create_window.php";
        </script>
        </body>
        </html>
<?php
	} else {
		//========================================================
		// Hapmap directory doesn't exist, go about creating it.
		//--------------------------------------------------------

		// Create the hapmap folder inside the user's hapmaps directory
		mkdir($hapmap_dir);
		secureNewDirectory($hapmap_dir);
		chmod($hapmap_dir,0777);

		// Initialize 'process_log.txt' file.
		$logOutputName = $hapmap_dir."/process_log.txt";
		$logOutput     = fopen($logOutputName, 'a');
		fwrite($logOutput, "Running 'scripts_seqModules/scripts_hapmaps/hapmap.install_3.php'.\n");
		fwrite($logOutput, "\tuser            = ".$user."\n");
		fwrite($logOutput, "\thapmap          = ".$hapmap."\n");
		fwrite($logOutput, "\tgenome          = ".$genome."\n");
		fwrite($logOutput, "\treferencePloidy = ".$referencePloidy."\n");
		fwrite($logOutput, "\tproject1        = ".$project1."\n");
		fwrite($logOutput, "\tproject2        = ".$project2."\n");
		fwrite($logOutput, "\tcolorA          = ".$colorA."\n");
		fwrite($logOutput, "\tcolorB          = ".$colorB."\n");

		// Create 'colors.txt' file to contain colors used in haplotype figures.
		$handleName = $hapmap_dir."/colors.txt";
		if (file_exists($handleName)) {
		} else {
			$handle     = fopen($handleName, 'w');
			fwrite($handle, $colorA."\n".$colorB);
			fclose($handle);
		}

		// Create 'genome.txt' file to contain genome used in haplotype.
		$handleName = $hapmap_dir."/genome.txt";
		if (file_exists($handleName)) {
		} else {
			$handle     = fopen($handleName, 'w');
			fwrite($handle, $genome);
			fclose($handle);
		}

		// Create 'parent.txt' file to contain parent genome used in haplotype.
		$handleName = $hapmap_dir."/parent.txt";
		if (file_exists($handleName)) {
		} else {
			$handle     = fopen($handleName, 'w');
			if ($referencePloidy == 2) {
				fwrite($handle, $project1);
			} else {
				fwrite($handle, $project1."\n".$project2);
			}
			fclose($handle);
		}

		// Initialize 'haplotype.txt' file to hold haplotype entry descriptions.
		if ($referencePloidy == 2) {
			$haplotypeFileName = $hapmap_dir."/haplotypeMap.txt";
			$haplotypeFile     = fopen($haplotypeFileName, 'a');
			fwrite($haplotypeFile, $project1."\n".$project2."\n".$hapmap_description."\n");
			fclose($haplotypeFile);
		} else {
			// for haplotype map derived from two haploid strains, there will not be more entries after construction.
			// Because of this, 'parent.txt' will contain sufficient information to define haplotype map.
		}

		// Generate 'working.txt' file to let pipeline know that processing is underway.
		$handleName      = $hapmap_dir."/working.txt";
		$handle          = fopen($handleName, 'w');
		$startTimeString = date("Y-m-d H:i:s");
		fwrite($handle, $startTimeString);
		fclose($handle);

		fwrite($logOutput, "'scripts_seqModules/scripts_hapmaps/hapmap.install_3.php' completed.\n");
		fwrite($logOutput, "...........................................................\n");
		fwrite($logOutput, $currentPath."\n");
		fclose($logOutput);

		// Pass control over to a shell script ('scripts_seqModules/scripts_hapmaps/hapmap.install_4.sh') to continue processing and link with matlab.
		$system_call_string = "sh hapmap.install_4.sh ".$user." ".$referencePloidy." ".$project1." ".$project2." ".$hapmap." > /dev/null &";
		system($system_call_string);
?>
<script type="text/javascript">
	var el3 = parent.document.getElementById('Hidden_GenerateNewHapmap');
	el3.style.display = 'none';

	var ff = parent.document.getElementById('panel_hapmap_iframe');
	ff.src = ff.src;
</script>
<?php
        }
?>
