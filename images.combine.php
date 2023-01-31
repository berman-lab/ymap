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
	$projectsShown = trim(filter_input(INPUT_POST, "projectsShown", FILTER_SANITIZE_STRING));	// strip out any html tags.
	$projectsShown = str_replace(" ","_",$projectsShown);						// convert any spaces to underlines.
	$projectsShown = preg_replace("/[\s\W]+/", "", $projectsShown);					// remove everything but alphanumeric characters and underlines.
	$projectsShown = str_replace("_"," ",$projectsShown);						// converts underlines back to spaces for later processing.


	// Script does not securely process received input from client.
	// It also doesn't produce intended output, so exiting now.
	exit();

	// auxillary functions
	// sets image background to white and init default image parameters
	function setBackgroundWhite($image) {
		$white = imagecolorallocate($image , 255, 255, 255);
		imagefill($image, 0, 0, $white);
		// setting default image parameters
		imagealphablending($image, false);
		imagesavealpha($image, true);
	}

	require_once 'constants.php';
	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";
	echo "projectsShown = '".$projectsShown."'<br><br>\n";
	// general variables
	$linearCartoonHeight = 140; // the height in px of the cartoon without labels 136 valid so 4px
	// creating combined images
	$projectsShown_entries = explode(" ",$projectsShown);
	$initial_entry = $projectsShown_entries[0];
	$entry_parts   = explode(":",$initial_entry);
	$fig_user      = $user;
	$fig_project   = $entry_parts[1];
	$fig_key       = $entry_parts[2];
	$fig_color1    = $entry_parts[3];
	$fig_color2    = $entry_parts[4];
	$fig_parent    = $entry_parts[5];
	$fig_CNV_SNP   = "users/".$fig_user."/projects/".$fig_project."/fig.CNV-SNP-map.2.png";
	$fig_CNV       = "users/".$fig_user."/projects/".$fig_project."/fig.CNV-map.2.png";
	if (file_exists("users/".$fig_user."/projects/".$fig_project."/fig.allelic_ratio-map.c2.png")) {   // ddRADseq.
		$fig_SNP   = "users/".$fig_user."/projects/".$fig_project."/fig.allelic_ratio-map.c2.png";
	} else { // other.
		$fig_SNP   = "users/".$fig_user."/projects/".$fig_project."/fig.SNP-map.2.png";
	}
	if (file_exists($fig_CNV_SNP)) {     $initial_image = $fig_CNV_SNP;
	} elseif (file_exists($fig_CNV)) {   $initial_image = $fig_CNV;
	} elseif (file_exists($fig_SNP)) {   $initial_image = $fig_SNP;
	}
	$image_size    = getimagesize($initial_image);
	$image_width   = $image_size[0];
	$image_height  = $image_size[1];
	$numImages     = count($projectsShown_entries);
	// creating new images containers, width remains the same, height the same for first picture and then 130px for every other image (to contain only the cartoons)
	// and setting backgorund to white
	$working1      = imagecreatetruecolor($image_width,$image_height + ($numImages - 1)*$linearCartoonHeight);
	setBackgroundWhite($working1); // set background to white
	$working2      = imagecreatetruecolor($image_width,$image_height + ($numImages - 1)*$linearCartoonHeight);
	setBackgroundWhite($working2); // set background to white
	$working3      = imagecreatetruecolor($image_width,$image_height + ($numImages - 1)*$linearCartoonHeight);
	setBackgroundWhite($working3); // set background to white
	// connecting pictures
	foreach ($projectsShown_entries as $entry_key => $projectsShown_entry) {
		$entry_parts = explode(":",$projectsShown_entry);
		$fig_user    = $entry_parts[0];
		$fig_project = $entry_parts[1];
		$fig_key     = $entry_parts[2];
		$fig_color1  = $entry_parts[3];
		$fig_color2  = $entry_parts[4];
		$fig_parent  = $entry_parts[5];
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
		$image1_height  = $image_size[1];
		$image2_size    = getimagesize($fig_CNV);
		$image2_height  = $image_size[1];
		$image3_size    = getimagesize($fig_SNP);
		$image3_height  = $image_size[1];
		// loading white color for in between filling
		$white = imagecolorallocate($working1, 255, 255, 255);
		if ($entry_key == 0) {
			// copy first picture entirely
			imagecopy($working1, $image1, 0, 0, 0, 0,  $image_width, $image_height);
			imagecopy($working2, $image2, 0, 0, 0, 0,  $image_width, $image_height);
			imagecopy($working3, $image3, 0, 0, 0, 0,  $image_width, $image_height);
			echo "[] ".$entry_key." ".$fig_CNV_SNP."<br>\n";
		} else {
			// copy only cartoon
			imagecopy($working1, $image1, 0, $image_height + $linearCartoonHeight*($entry_key-1), 0, $image1_height - $linearCartoonHeight, $image_width, $linearCartoonHeight);
			// filing in white in between
			imagefilledrectangle($working1, 50, $image_height + $linearCartoonHeight*($entry_key-1) - 10, $image_width,  $image_height + $linearCartoonHeight*($entry_key-1) + 4, $white);
			imagecopy($working2, $image2, 0, $image_height + $linearCartoonHeight*($entry_key-1), 0, $image2_height - $linearCartoonHeight, $image_width, $linearCartoonHeight);
			// filing in white in between
			imagefilledrectangle($working2, 50, $image_height + $linearCartoonHeight*($entry_key-1)-10, $image_width,  $image_height + $linearCartoonHeight*($entry_key-1) + 4, $white);		
			imagecopy($working3, $image3, 0, $image_height + $linearCartoonHeight*($entry_key-1), 0, $image3_height - $linearCartoonHeight, $image_width, $linearCartoonHeight);
			// filing in white in between
			imagefilledrectangle($working3, 50, $image_height + $linearCartoonHeight*($entry_key-1) - 10, $image_width,  $image_height + $linearCartoonHeight*($entry_key-1) + 4, $white);			
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
