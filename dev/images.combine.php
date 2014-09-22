<?php
	session_start();
	if(!isset($_SESSION['logged_on'])){ ?> <script type="text/javascript"> parent.reload(); </script><?php echo "\n"; } else { $user = $_SESSION['user']; }
	require_once 'php/constants.php';
	echo "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01 Transitional//EN\" \"http://www.w3.org/TR/html4/loose.dtd\">\n";

	$projectsShown = filter_input(INPUT_POST, "projectsShown", FILTER_SANITIZE_STRING);
//	$user          = 'darren1';
//	$projectsShown = 'darren1:12353_A21-s02-m08-r09:1:null:null:SC5314_A21-s02-m08-r09 darren1:12353_vs_hapmap_A21-s02-m08-r09:3:cyan:magenta:12353_vs_hapmap_A21-s02-m08-r09 darren1:testing_hapmap:18:green:blue:testing_hapmap darren1:aaaaa:14:null:null:none default:Fig_02A.YQ2_raw:19:null:null:Fig_02A darren1:Fig_04_11461:5:null:null:Fig_04_11461';

	echo "projectsShown = '".$projectsShown."'<br><br>\n";

	$projectsShown_entries = explode(" ",$projectsShown);
	$initial_entry = $projectsShown_entries[0];
	$entry_parts   = explode(":",$initial_entry);
	$fig_user      = $entry_parts[0];
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
	$image_width   = $image_size[0];           // 1152
	$image_height  = $image_size[1];           // 90
	$numImages     = count($projectsShown_entries);
	$working1      = imagecreatetruecolor($image_width,$image_height*$numImages - 25*($numImages-1) - 5*$numImages);
	$working2      = imagecreatetruecolor($image_width,$image_height*$numImages - 25*($numImages-1) - 5*$numImages);
	$working3      = imagecreatetruecolor($image_width,$image_height*$numImages - 25*($numImages-1) - 5*$numImages);
	imagealphablending($working1, false);
	imagealphablending($working2, false);
	imagealphablending($working3, false);
	imagesavealpha(    $working1, true);
	imagesavealpha(    $working2, true);
	imagesavealpha(    $working3, true);
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
		$image1      = imagecreatefrompng($fig_CNV_SNP);
		$image2      = imagecreatefrompng($fig_CNV);
		$image3      = imagecreatefrompng($fig_SNP);
		if ($entry_key == 0) {
			imagecopy(    $working1, $image1, 0, 0,              0, 0,  $image_width, $image_height);
			imagecopy(    $working2, $image2, 0, 0,              0, 0,  $image_width, $image_height);
			imagecopy(    $working3, $image3, 0, 0,              0, 0,  $image_width, $image_height);
			echo "[] ".$entry_key." ".$fig_CNV_SNP."<br>\n";
		} else {
			imagecopy(    $working1, $image1, 0, 85+60*($entry_key-1), 0, 25, $image_width, $image_height);
			imagecopy(    $working2, $image2, 0, 85+60*($entry_key-1), 0, 25, $image_width, $image_height);
			imagecopy(    $working3, $image3, 0, 85+60*($entry_key-1), 0, 25, $image_width, $image_height);
			echo "...   ".$entry_key." ".$fig_CNV_SNP."<br>\n";
		}
	}
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
