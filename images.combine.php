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
	$projectsShown = sanitizeProjectsShown_POST("projectsShown");

	// auxillary functions
	// sets image background to white and init default image parameters
	function setBackgroundWhite($image) {
		$white = imagecolorallocate($image , 255, 255, 255);
		imagefill($image, 0, 0, $white);
		// setting default image parameters
		imagealphablending($image, false);
		imagesavealpha($image, true);
	}

	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";
	echo "projectsShown = '".$projectsShown."'<br><br>\n";

	// general variables
	$linearCartoonHeight = 139; //139 the height in px of the cartoon without labels 136 valid so + 4px.

	// break projectsShown string into individual strings per project.
	$projectsShown_entries = explode(" ",$projectsShown);

	// process first project string;
	$initial_entry = $projectsShown_entries[0];
	$entry_parts   = explode(":",$initial_entry);
	$fig_user      = $entry_parts[0];
	$fig_project   = $entry_parts[1];
	$fig_key       = $entry_parts[2];

	//================================================================
	// Validate sub-strings received from first projectsShown string.
	//----------------------------------------------------------------
	// Confirm requested user exists.
	$user_dir = "users/".$fig_user;
	if (!is_dir($user_dir)) {
		// user doesn't exist, should never happen: Force logout.
		session_destroy();
		?><script type="text/javascript"> parent.location.reload(); </script><?php
	}
	// Confirm requested project exists.
	$project_dir = "users/".$fig_user."/projects/".$fig_project;
	if (!is_dir($project_dir)) {
		// project doesn't exist, should never happen: Force logout.
		session_destroy();
		?><script type="text/javascript"> parent.location.reload(); </script><?php
	}

	// Determine initial figure strings.
	$fig_CNV_SNP   = "users/".$fig_user."/projects/".$fig_project."/fig.CNV-SNP-map.2.png";
	$fig_CNV       = "users/".$fig_user."/projects/".$fig_project."/fig.CNV-map.2.png";
	if (file_exists("users/".$fig_user."/projects/".$fig_project."/fig.allelic_ratio-map.c2.png")) {   // ddRADseq.
		$fig_SNP   = "users/".$fig_user."/projects/".$fig_project."/fig.allelic_ratio-map.c2.png";
	} else { // other.
		$fig_SNP   = "users/".$fig_user."/projects/".$fig_project."/fig.SNP-map.2.png";
	}
	if (file_exists($fig_CNV_SNP)) { $initial_image = $fig_CNV_SNP; }
	elseif (file_exists($fig_CNV)) { $initial_image = $fig_CNV; }
	elseif (file_exists($fig_SNP)) { $initial_image = $fig_SNP; }
	$image_size    = getimagesize($initial_image);
	$image_width   = $image_size[0];
	$image_height  = $image_size[1];

	// determine number of images to combine.
	$numImages     = count($projectsShown_entries);

	// Make new images containers; width is same for all; height is same for first, then +140px for others (to contain only the cartoons).
	$working1      = imagecreatetruecolor($image_width,$image_height + ($numImages - 1)*$linearCartoonHeight);
	$working2      = imagecreatetruecolor($image_width,$image_height + ($numImages - 1)*$linearCartoonHeight);
	$working3      = imagecreatetruecolor($image_width,$image_height + ($numImages - 1)*$linearCartoonHeight);

	// Set backgrounds to white.
	setBackgroundWhite($working1);
	setBackgroundWhite($working2);
	setBackgroundWhite($working3);

	//=======================================================
	// Clean up any previously constructed combined figures.
	//-------------------------------------------------------
	$image_file1 = "users/".$user."/combined_figure.1.png";
	$image_file2 = "users/".$user."/combined_figure.2.png";
	$image_file3 = "users/".$user."/combined_figure.3.png";
	if (file_exists($image_file1)) { unlink($image_file1); }
	if (file_exists($image_file2)) { unlink($image_file2); }
	if (file_exists($image_file3)) { unlink($image_file3); }

	//=========================
	// Build combined figures.
	//-------------------------
	foreach ($projectsShown_entries as $entry_key => $projectsShown_entry) {
		// process subsequent projectsShown strings.
		$entry_parts = explode(":",$projectsShown_entry);
		$fig_user    = $entry_parts[0];
		$fig_project = $entry_parts[1];
		$fig_key     = $entry_parts[2];

		//================================================================
		// Validate sub-strings received from projectsShown string.
		//----------------------------------------------------------------
		// Confirm requested user exists.
		$user_dir = "users/".$fig_user;
		if (!is_dir($user_dir)) {
			// user doesn't exist, should never happen: Force logout.
			session_destroy();
			?><script type="text/javascript"> parent.location.reload(); </script><?php
		}
		// Confirm requested project exists.
		$project_dir = "users/".$fig_user."/projects/".$fig_project;
		if (!is_dir($project_dir)) {
			// project doesn't exist, should never happen: Force logout.
			session_destroy();
			?><script type="text/javascript"> parent.location.reload(); </script><?php
		}

		// Determine figure strings.
		$fig_CNV_SNP = "users/".$fig_user."/projects/".$fig_project."/fig.CNV-SNP-map.2.png";
		$fig_CNV     = "users/".$fig_user."/projects/".$fig_project."/fig.CNV-map.2.png";
		if (file_exists("users/".$fig_user."/projects/".$fig_project."/fig.allelic_ratio-map.c2.png")) {   // ddRADseq.
			$fig_SNP = "users/".$fig_user."/projects/".$fig_project."/fig.allelic_ratio-map.c2.png";
		} else { // other.
			$fig_SNP = "users/".$fig_user."/projects/".$fig_project."/fig.SNP-map.2.png";
		}

		// creating images for copy
		$image1      = imagecreatefrompng($fig_CNV_SNP);
		$image2      = imagecreatefrompng($fig_CNV);
		$image3      = imagecreatefrompng($fig_SNP);

		// getting sizes
		$image1_size    = getimagesize($fig_CNV_SNP);
		$image2_size    = getimagesize($fig_CNV);
		$image3_size    = getimagesize($fig_SNP);
		$image1_height  = $image_size[1];
		$image2_height  = $image_size[1];
		$image3_height  = $image_size[1];

		// loading white color for in between filling
		$white = imagecolorallocate($working1, 255, 255, 255);

		if ($entry_key == 0) {
			// copy first figure entirely
			imagecopy($working1, $image1, 0, 0, 0, 0,  $image_width, $image_height);
			imagecopy($working2, $image2, 0, 0, 0, 0,  $image_width, $image_height);
			imagecopy($working3, $image3, 0, 0, 0, 0,  $image_width, $image_height);

			echo "[] ".$entry_key." ".$fig_CNV_SNP."<br>\n";
		} else {
			// copy individual images together.
			imagecopy($working1, $image1, 0, $image_height + $linearCartoonHeight*($entry_key-1), 0, $image1_height - $linearCartoonHeight, $image_width, $linearCartoonHeight);
			imagecopy($working2, $image2, 0, $image_height + $linearCartoonHeight*($entry_key-1), 0, $image2_height - $linearCartoonHeight, $image_width, $linearCartoonHeight);
			imagecopy($working3, $image3, 0, $image_height + $linearCartoonHeight*($entry_key-1), 0, $image3_height - $linearCartoonHeight, $image_width, $linearCartoonHeight);

			// fill white in between.
			imagefilledrectangle($working1, 50, $image_height + $linearCartoonHeight*($entry_key-1)-10, $image_width,  $image_height + $linearCartoonHeight*($entry_key-1)+4, $white);
			imagefilledrectangle($working2, 50, $image_height + $linearCartoonHeight*($entry_key-1)-10, $image_width,  $image_height + $linearCartoonHeight*($entry_key-1)+4, $white);
			imagefilledrectangle($working3, 50, $image_height + $linearCartoonHeight*($entry_key-1)-10, $image_width,  $image_height + $linearCartoonHeight*($entry_key-1)+4, $white);

			echo "...   ".$entry_key." ".$fig_CNV_SNP."<br>\n";
		}
	}

// saving images and destroying
imagepng($working1,"users/".$user."/combined_figure.1.png");
imagepng($working2,"users/".$user."/combined_figure.2.png");
imagepng($working3,"users/".$user."/combined_figure.3.png");
imagedestroy($working1);
imagedestroy($working2);
imagedestroy($working3);
imagedestroy($image1);
imagedestroy($image2);
imagedestroy($image3);
?>
