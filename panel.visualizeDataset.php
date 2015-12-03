<?php
	session_start();
	if (isset($_SESSION['logged_on'])) { $user = $_SESSION['user']; } else { $user = 'default'; }
?>
<style type="text/css">
	html * {
		font-family: arial !important;
	}
</style>
<font size='3'>View figures for installed datasets at bottom of page by selecting checkboxes.</font><br><br>
<table width="100%" cellpadding="0"><tr>
<td width="65%" valign="top">
	<?php
	//.---------------.
	//| User projects |
	//'---------------'
	if (isset($_SESSION['logged_on'])) {
		$projectsDir      = "users/".$user."/projects/";
		$projectFolders   = array_diff(glob($projectsDir."*"), array('..', '.'));
		// Sort directories by date, newest first.
		array_multisort(array_map('filemtime', $projectFolders), SORT_DESC, $projectFolders);
		// Trim path from each folder string.
		foreach($projectFolders as $key=>$folder) {   $projectFolders[$key] = str_replace($projectsDir,"",$folder);   }
		// Split project list into ready/working/starting lists for sequential display.
		$projectFolders_complete = array();
		$projectFolders_working  = array();
		$projectFolders_starting = array();
		foreach($projectFolders as $key=>$project) {
			if (file_exists("users/".$user."/projects/".$project."/complete.txt")) {
				array_push($projectFolders_complete,$project);
			} else if (file_exists("users/".$user."/projects/".$project."/working.txt")) {
				array_push($projectFolders_working, $project);
			} else {
				array_push($projectFolders_starting,$project);
			}
		}
		$userProjectCount_starting = count($projectFolders_starting);
		$userProjectCount_working  = count($projectFolders_working);
		$userProjectCount_complete = count($projectFolders_complete);

		// Sort complete and working projects alphabetically.
		array_multisort($projectFolders_working,  SORT_ASC, $projectFolders_working);
		array_multisort($projectFolders_complete, SORT_ASC, $projectFolders_complete);
		// Build new 'projectFolders' array;
		$projectFolders   = array();
		$projectFolders   = array_merge($projectFolders_starting, $projectFolders_working, $projectFolders_complete);
		$userProjectCount = count($projectFolders);

		echo "<b><font size='2'>User installed datasets:</font></b>\n\t\t";
		echo "<br>\n\t\t";
		foreach($projectFolders_starting as $key_=>$project) {
			// Load colors for project.
			$colors_file  = "users/".$user."/projects/".$project."/colors.txt";
			if (file_exists($colors_file)) {
				$handle       = fopen($colors_file,'r');
				$colorString1 = trim(fgets($handle));
				$colorString2 = trim(fgets($handle));
				fclose($handle);
			} else {
				$colorString1 = 'null';
				$colorString2 = 'null';
			}

			$dataType_file        = "users/".$user."/projects/".$project."/dataType.txt";
			if (file_exists($dataType_file)) {
				$handle       = fopen($dataType_file,'r');
				$dataType     = trim(fgets($handle));
				fclose($handle);
			} else {
				$dataType     = 'null';
			}
			if (strcmp($dataType,"0") == 0) {
				$colorString1 = "cyan";
				$colorString2 = "magenta";
			}

			$parent_file          = "users/".$user."/projects/".$project."/parent.txt";
			$handle               = fopen($parent_file,'r');
			$parentString         = trim(fgets($handle));
			fclose($handle);
			$key = $key_;
			echo "<span id='p_label_".$key."' style='color:#CC0000;'>\n\t\t";
			echo "<font size='2'>".($key+1).".";
			echo "<input id='show_".$key."' type='checkbox' onclick=\"parent.openProject('".$user."','".$project."','".$key."','".$colorString1."','".$colorString2."','".$parentString."');\" style=\"visibility:hidden;\">";
			echo "\n\t\t".$project."</font></span>\n\t\t";
			echo "<span id='p_".$project."_type'></span>\n\t\t";
			echo "<br>\n\t\t";
			echo "<div id='frameContainer.p2_".$key."'></div>";
		}
		foreach($projectFolders_working as $key_=>$project) {
			// Load colors for project.
			$colors_file          = "users/".$user."/projects/".$project."/colors.txt";
			if (file_exists($colors_file)) {
				$handle       = fopen($colors_file,'r');
				$colorString1 = trim(fgets($handle));
				$colorString2 = trim(fgets($handle));
				fclose($handle);
			} else {
				$colorString1 = 'null';
				$colorString2 = 'null';
			}

			$dataType_file        = "users/".$user."/projects/".$project."/dataType.txt";
			if (file_exists($dataType_file)) {
				$handle       = fopen($dataType_file,'r');
				$dataType     = trim(fgets($handle));
				fclose($handle);
			} else {
				$dataType     = 'null';
			}
			if (strcmp($dataType,"0") == 0) {
				$colorString1 = "cyan";
				$colorString2 = "magenta";
			}

			$parent_file          = "users/".$user."/projects/".$project."/parent.txt";
			$handle               = fopen($parent_file,'r');
			$parentString         = trim(fgets($handle));
			fclose($handle);
			$key = $key_ + $userProjectCount_starting;
			echo "<span id='p_label_".$key."' style='color:#BB9900;'>\n\t\t";
			echo "<font size='2'>".($key+1).".";
			echo "<input  id='show_".$key."' type='checkbox' onclick=\"parent.openProject('".$user."','".$project."','".$key."','".$colorString1."','".$colorString2."','".$parentString."');\" style=\"visibility:hidden;\">";
			echo "\n\t\t".$project."</font></span>\n\t\t";
			echo "<span id='p_".$project."_type'></span>\n\t\t";
			echo "<br>\n\t\t";
			echo "<div id='frameContainer.p2_".$key."'></div>";
		}
		foreach($projectFolders_complete as $key_=>$project) {
			// Load colors for project.
			$colors_file          = "users/".$user."/projects/".$project."/colors.txt";
			if (file_exists($colors_file)) {
				$handle       = fopen($colors_file,'r');
				$colorString1 = trim(fgets($handle));
				$colorString2 = trim(fgets($handle));
				fclose($handle);
			} else {
				$colorString1 = 'null';
				$colorString2 = 'null';
			}

			$dataType_file        = "users/".$user."/projects/".$project."/dataType.txt";
			if (file_exists($dataType_file)) {
				$handle       = fopen($dataType_file,'r');
				$dataType     = trim(fgets($handle));
				fclose($handle);
			} else {
				$dataType     = 'null';
			}
			if (strcmp($dataType,"0") == 0) {
				$colorString1 = "cyan";
				$colorString2 = "magenta";
			}
			
			$json_file_list = json_encode(scandir("users/$user/projects/$project"));
			$parent_file          = "users/".$user."/projects/".$project."/parent.txt";
			$handle               = fopen($parent_file,'r');
			$parentString         = trim(fgets($handle));
			fclose($handle);
			$key = $key_ + $userProjectCount_starting + $userProjectCount_working;
			echo "<span id='project_label_".$key."' style='color:#00AA00;'>\n\t\t";
			echo "<font size='2'>".($key+1).".";
			echo "<input  id='show_$key' type='checkbox' onclick=\"parent.openProject('$user','$project','$key','$colorString1','$colorString2','$parentString');\" data-file-list='$json_file_list' >";
			echo "\n\t\t".$project."</font></span>\n\t\t";
			echo "<span id='p2_".$project."_delete'></span><span id='p_".$project."_type'></span>\n\t\t";
			echo "<br>\n\t\t";
			echo "<div id='frameContainer.p1_".$key."'></div>";
		}
	} else {
		$userProjectCount_starting = 0;
		$userProjectCount_working  = 0;
		$userProjectCount_complete = 0;
	}
	?>
</td><td width="35%" valign="top">
	<br><?php
	//.-----------------.
	//| System projects |
	//'-----------------'
	$projectsDir          = "users/default/projects/";
	$systemProjectFolders = array_diff(glob($projectsDir."*"), array('..', '.'));
	// Sort directories by date, newest first.
	array_multisort($systemProjectFolders, SORT_ASC, $systemProjectFolders);
	// array_multisort(array_map('filemtime', $systemProjectFolders), SORT_DESC, $systemProjectFolders);
	// Trim path from each folder string.
	foreach($systemProjectFolders as $key=>$folder) {   $systemProjectFolders[$key] = str_replace($projectsDir,"",$folder);   }
	$systemProjectCount = count($systemProjectFolders);
	echo "<b><font size='2'>Sample datasets:</font></b>\n\t\t";
	echo "<br>\n\t\t";
	foreach ($systemProjectFolders as $key_=>$project) {
		// Load colors for project.
		$colors_file          = "users/default/projects/".$project."/colors.txt";
		if (file_exists($colors_file)) {
			$handle       = fopen($colors_file,'r');
			$colorString1 = trim(fgets($handle));
			$colorString2 = trim(fgets($handle));
			fclose($handle);
		} else {
			$colorString1 = 'null';
			$colorString2 = 'null';
		}
		$json_file_list = json_encode(scandir("users/default/projects/$project"));
		$parent_file  = "users/default/projects/".$project."/parent.txt";
		$handle       = fopen($parent_file,'r');
		$parentString = trim(fgets($handle));
		fclose($handle);

		$key = $key_ + $userProjectCount_starting + $userProjectCount_working + $userProjectCount_complete;
		echo "<font size='2'>".($key+1).".";
		echo "<input  id='show_".$key."_sys' type='checkbox' onclick=\"parent.openProject('default','".$project."','".$key."_sys','".$colorString1."','".$colorString2."','".$parentString."');\" data-file-list='".$json_file_list. "'>";
		echo $project."</font>";
		echo "<br>\n\t\t";
	}
	?>
</td></tr></table>
<script type="text/javascript">

if(localStorage.getItem("projectsShown")){
	var projectsShown = localStorage.getItem("projectsShown");
}

function Display_sample_figures() {
	localStorage.setItem("projectsShown","");
<?php
	if (isset($_SESSION['logged_on'])) {
		foreach ($systemProjectFolders as $key=>$project) {
			if ($key < 2) {
				$new_key = $key+$userProjectCount; // offset example datasets by number of user projects.
				echo "\tvar show_button_element = document.getElementById('show_".$key."_sys');\n";
				echo "\tshow_button_element.checked = true;\n";
				echo "\tparent.openProject('default','".$project."','".$new_key."_sys','null','null','null');\n";
			}
		}
	}
	?>
}
Display_sample_figures();

</script>
