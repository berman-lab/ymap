<?php
	session_start();
	require_once 'php/constants.php';
?>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" "http://www.w3.org/TR/html4/loose.dtd">
<html lang="en">
	<head>
		<style type="text/css">
			body {
				font-family: arial;
			}
			.upload {
				width:          100%;
				border:         0;
				height:         0px;
				vertical-align: middle;
				align:          left;
			}
			.tab {
				margin-left:    1cm;
			}
		</style>
		<meta http-equiv="content-type" content="text/html; charset=utf-8">
		<title>Y-MAP</title>
<!-- Used by secondary pages to update page on completion of processing. --!>
		<script type="text/javascript">
		resize_iframe = function(iframe_key, pixels) {
			document.getElementById(iframe_key).style.height = pixels+"px";
		}
		update_project_remove_iframe = function(project_key) {
			project_key                          = project_key.replace('p_','');
			var show_button_element              = document.getElementById('show_'+project_key);
            show_button_element.style.visibility = 'visible';
			var project_iframe                   = document.getElementById('frameContainer.p2_'+project_key);
			project_iframe.innerHTML             = '';
			document.getElementById("p_".project_key).style.height = "0px";
		}
		update_project_label_color = function(project_key,label_color) {
			project_key               = project_key.replace('p_','');
			var project_label         = document.getElementById('project_label_'+project_key);
			project_label.style.color = label_color;
		}
		update_hapmap_label_color = function(hapmap_key,label_color) {
			hapmap_key               = hapmap_key.replace('h_','');
			var hapmap_label         = document.getElementById('hapmap_label_'+hapmap_key);
			hapmap_label.style.color = label_color;
		}
		update_genome_label_color = function(genome_key,label_color) {
			genome_key               = genome_key.replace('g_','');
			var genome_label         = document.getElementById('genome_label_'+genome_key);
			genome_label.style.color = label_color;
		}
		resize_project = function(project_key, pixels) {
			document.getElementById("p_".project_key).style.height = pixels+"px";
		}
		resize_hapmap  = function(hapmap_key, pixels) {
			document.getElementById("h_".project_key).style.height = pixels+"px";
		}
		resize_genome  = function(genome_key, pixels) {
			document.getElementById("g_".project_key).style.height = pixels+"px";
		}
		var user = "<?php echo $_SESSION['user']?>";
		</script>

	</head>
	<body>
<table width="100%"><tr>
	<td width="25%" align="center" style="max-height:100%">
		<img src="images/Logo_title.2.png" alt="Y-MAP; Yeast Mapping Analysis Pipeline">
	</td><td width="75%" align="right" valign="top">

<!---------------------------------------------------------------------------!>
<!--- Construct tabbed menu and content pane.  ------------------------------!>
<!---------------------------------------------------------------------------!>
		<table width="100%" height="210px" cellspacing="0">
		<tr>
			<td valign="bottom" style="height:15px; width:80px;" align="center" id="tab_user"     ><div onclick="tabWindow_user();"     >User</div></td>
			<td valign="bottom" style="height:15px; width:80px;" align="center" id="tab_project"  ><div onclick="tabWindow_project();"  >Project</div></td>
			<td valign="bottom" style="height:15px; width:80px;" align="center" id="tab_genome"   ><div onclick="tabWindow_genome();"   >Genome</div></td>
			<td valign="bottom" style="height:15px; width:80px;" align="center" id="tab_hapmap"   ><div onclick="tabWindow_hapmap();"   >Hapmap</div></td>
			<td valign="bottom" style="height:15px; width:80px;" align="center" id="tab_system"   ><div onclick="tabWindow_system();"   >System</div></td>
			<td valign="bottom" style="height:15px; width:80px;" align="center" id="tab_about"    ><div onclick="tabWindow_about();"    >About</div></td>
			<td valign="bottom" style="height:15px; width:160px;" align="center" id="tab_examples" ><div onclick="tabWindow_examples();" >Example Datasets</div></td>
			<td valign="bottom" id="tab_blank">&nbsp;</td>
		</tr><tr>
			<td colspan="8" valign="top" id="tab_content">
<!---------------------------------------------------------------------------!>
			<div id="panel_user" name="panel_user">
				<?php
				if (isset($_SESSION['logged_on'])) {
					echo "User '<b>".$_SESSION['user']."</b>' logged in. ";
					// provide logout button.
					$user = $_SESSION['user'];
					echo "<button type='button' onclick=\"update_projectsShown_after_logout(); window.location.href='php/logout_server.php'\">Logout</button>\n\t\t\t\t";
					// provide delete-user button.
					echo "<span id='u_".$user."_delete'></span>\n\t\t\t\t";
					echo "<span id='deleteUser'><button type='button' onclick=\"deleteUserConfirmation('".$user."','deleteUser')\">Delete User</button></span>\n\t\t\t\t";
					echo "<br><br>\n\t\t\t";
					echo "<font size='2'>";
					echo "You can navigate through the above menu and show/close projects while new datafiles are uploading.<br>";
					echo "A page reload or project/genome/hapmap creation/deletion, however, will interrupt file transfer.<br><br>";
					echo "Depending on system load, tasks may take an hour or more to complete after data upload is complete.<br><br>";
					echo "Reload page and select 'projects' tab to check for newly completed projects.";
					echo "</font>";
				} else {
					echo "<script type=\"text/javascript\">update_projectsShown_after_logout();</script>";
					echo "<form action='php/login_server.php' method='post'>";
					echo "<label for='user'>Username: </label><input type='text' id='user' name='user'><br>";
					echo "<label for='pw'>Password: </label><input type='password' id='pw' name='pw'><br>";
					echo "<button type='submit' onclick=\"update_projectsShown_after_logout();\">Login</button>";
					echo "</form><br>";
					echo "<font size='2'>";
					echo "If you don't have a user account, you may make one by clicking below.<br>";
					echo "<button type='button' onclick=\"update_projectsShown_after_logout(); window.location.href='register.php'\">Register new user.</button>";
					echo "</font>";
				}
			?></div>
<!---------------------------------------------------------------------------!>
			<div id="panel_project" name="panel_project">
				<table width="100%" cellpadding="0"><tr><td width="65%">
				<?php
				//.---------------.
                //| User projects |
                //'---------------'
				$userProjectCount = 0;
				if (isset($_SESSION['logged_on'])) {
					$projectsDir    = $directory."users/".$user."/projects/";
					$projectFolders = array_diff(glob($projectsDir."*"), array('..', '.'));
					// Sort directories by date, newest first.
					array_multisort(array_map('filemtime', $projectFolders), SORT_DESC, $projectFolders);
					// Trim path from each folder string.
					foreach($projectFolders as $key=>$folder) {   $projectFolders[$key] = str_replace($projectsDir,"",$folder);   }

					// Split project list into ready/working/starting lists for sequential display.
					$projectFolders_complete = array();
					$projectFolders_working  = array();
					$projectFolders_starting = array();
					foreach($projectFolders as $key=>$project) {
						if (file_exists($directory."users/".$user."/projects/".$project."/complete.txt")) {
							array_push($projectFolders_complete,$project);
						} else if (file_exists($directory."users/".$user."/projects/".$project."/working.txt")) {
							array_push($projectFolders_working, $project);
						} else {
							array_push($projectFolders_starting,$project);
						}
					}
					$userProjectCount_starting = count($projectFolders_starting);
					$userProjectCount_working  = count($projectFolders_working);
					$userProjectCount_complete = count($projectFolders_complete);
					$userProjectCount          = count($projectFolders);

					// Sort complete and working projects alphabetically.
					array_multisort($projectFolders_working,  SORT_ASC, $projectFolders_working);
					array_multisort($projectFolders_complete, SORT_ASC, $projectFolders_complete);

					// Build new 'projectFolders' array;
					$projectFolders   = array();
					$projectFolders   = array_merge($projectFolders_starting, $projectFolders_working, $projectFolders_complete);
					$userProjectCount = count($projectFolders);

					echo "<button type='button' onclick=\"update_projectsShown_after_new_project(); window.location.href='project.create.php';\">Install New Project</button><br>\n\t\t\t\t";
					echo "<b><font size='2'>User installed projects:</font></b>\n\t\t\t\t";
					echo "<br>\n\t\t\t\t";
					echo "<div class='tab' style='height:141px; overflow-y:scroll;'>\n\t\t\t\t";
					foreach($projectFolders_starting as $key_=>$project) {
						// Load colors for project.
						$colorFile        = $GLOBALS['directory']."users/".$user."/projects/".$project."/colors.txt";
						if (file_exists($colorFile)) {
							$handle       = fopen($colorFile,'r');
							$colorString1 = trim(fgets($handle));
							$colorString2 = trim(fgets($handle));
							fclose($handle);
						} else {
							$colorString1 = 'null';
							$colorString2 = 'null';
						}
						$parentFile   = $GLOBALS['directory']."users/".$user."/projects/".$project."/parent.txt";
						$handle       = fopen($parentFile,'r');
						$parentString = trim(fgets($handle));
						fclose($handle);

						$key = $key_;
						echo "<span id='project_label_".$key."' style='color:#CC0000;'>\n\t\t\t\t";
						echo "<font size='2'>".($key+1).". <button id='show_".$key."' type='button' onclick=\"openProject('".$user."','".$project."','".$key."','".$colorString1."','".$colorString2."','".$parentString."')\">Show</button>\n\t\t\t\t";
						echo $project."</font></span>\n\t\t\t\t";

						echo "<span id='p1_".$key."_delete'><button id='delete_".$key."' type='button' onclick=\"deleteProjectConfirmation('".$project."','p1_".$key."_delete')\">Delete</button>";
						echo "</span>\n\t\t\t\t";
						echo "<span id='p2_".$project."_delete'></span><span id='p_".$project."_type'></span>\n\t\t\t\t";
						echo "<br>\n\t\t\t\t";
						echo "<div id='frameContainer.p3_".$key."'></div>";
					}
					foreach($projectFolders_working as $key_=>$project) {
						// Load colors for project.
						$colorFile        = $GLOBALS['directory']."users/".$user."/projects/".$project."/colors.txt";
						if (file_exists($colorFile)) {
							$handle       = fopen($colorFile,'r');
							$colorString1 = trim(fgets($handle));
							$colorString2 = trim(fgets($handle));
							fclose($handle);
						} else {
							$colorString1 = 'null';
							$colorString2 = 'null';
						}
						$parentFile   = $GLOBALS['directory']."users/".$user."/projects/".$project."/parent.txt";
						$handle       = fopen($parentFile,'r');
						$parentString = trim(fgets($handle));
						fclose($handle);

						$key = $key_ + $userProjectCount_starting;
						echo "<span id='project_label_".$key."' style='color:#BB9900;'>\n\t\t\t\t";
						echo "<font size='2'>".($key+1).". <button id='show_".$key."' type='button' onclick=\"openProject('".$user."','".$project."','".$key."','".$colorString1."','".$colorString2."','".$parentString."')\">Show</button>\n\t\t\t\t";
						echo $project."</font></span>\n\t\t\t\t";

						echo "<span id='p1_".$key."_delete'><button id='delete_".$key."' type='button' onclick=\"deleteProjectConfirmation('".$project."','p1_".$key."_delete')\">Delete</button>";
						echo "</span>\n\t\t\t\t";
						echo "<span id='p2_".$project."_delete'></span><span id='p_".$project."_type'></span>\n\t\t\t\t";
						echo "<br>\n\t\t\t\t";
						echo "<div id='frameContainer.p2_".$key."'></div>";
					}
					foreach($projectFolders_complete as $key_=>$project) {
						// Load colors for project.
						$colorFile        = $GLOBALS['directory']."users/".$user."/projects/".$project."/colors.txt";
						if (file_exists($colorFile)) {
							$handle       = fopen($colorFile,'r');
							$colorString1 = trim(fgets($handle));
							$colorString2 = trim(fgets($handle));
							fclose($handle);
						} else {
							$colorString1 = 'null';
							$colorString2 = 'null';
						}
						$parentFile   = $GLOBALS['directory']."users/".$user."/projects/".$project."/parent.txt";
						$handle       = fopen($parentFile,'r');
						$parentString = trim(fgets($handle));
						fclose($handle);

						$key = $key_ + $userProjectCount_starting + $userProjectCount_working;
						echo "<span id='project_label_".$key."' style='color:#00AA00;'>\n\t\t\t\t";
						echo "<font size='2'>".($key+1).". <button id='show_".$key."' type='button' onclick=\"openProject('".$user."','".$project."','".$key."','".$colorString1."','".$colorString2."','".$parentString."')\">Show</button>\n\t\t\t\t";
						echo $project."</font></span>\n\t\t\t\t";

						echo "<span id='p1_".$key."_delete'><button id='delete_".$key."' type='button' onclick=\"deleteProjectConfirmation('".$project."','p1_".$key."_delete')\">Delete</button>";
						echo "</span>\n\t\t\t\t";
						echo "<span id='p2_".$project."_delete'></span><span id='p_".$project."_type'></span>\n\t\t\t\t";
						echo "<br>\n\t\t\t\t";
						echo "<div id='frameContainer.p1_".$key."'></div>";
					}
				}?></div></td><td width="35%"><br><?php
				//.-----------------.
				//| System projects |
				//'-----------------'
				$projectsDir          = $directory."users/default/projects/";
				$systemProjectFolders = array_diff(glob($projectsDir."*"), array('..', '.'));
				// Sort directories by date, newest first.
				array_multisort(array_map('filemtime', $systemProjectFolders), SORT_DESC, $systemProjectFolders);
				// Trim path from each folder string.
				foreach($systemProjectFolders as $key=>$folder) {   $systemProjectFolders[$key] = str_replace($projectsDir,"",$folder);   }
				$systemProjectCount = count($systemProjectFolders);
				echo "<b><font size='2'>System projects:</font></b>\n\t\t\t\t";
				echo "<br>\n\t\t\t\t";
				echo "<div class='tab' style='height:141px; overflow-y:scroll;'>\n\t\t\t\t";
				foreach ($systemProjectFolders as $key_=>$project) {
					// Load colors for project.
					$colorFile        = $GLOBALS['directory']."users/default/projects/".$project."/colors.txt";
					if (file_exists($colorFile)) {
						$handle       = fopen($colorFile,'r');
						$colorString1 = trim(fgets($handle));
						$colorString2 = trim(fgets($handle));
						fclose($handle);
					} else {
						$colorString1 = 'null';
						$colorString2 = 'null';
					}
					$parentFile   = $GLOBALS['directory']."users/default/projects/".$project."/parent.txt";
					$handle       = fopen($parentFile,'r');
					$parentString = trim(fgets($handle));
					fclose($handle);

					$key = $key_ + $userProjectCount_starting + $userProjectCount_working + $userProjectCount_complete;
					echo "<font size='2'>".($key+1).". <button id='show_".$key."' type='button' onclick=\"openProject('default','".$project."','".$key."','".$colorString1."','".$colorString2."','".$parentString."')\">Show</button>\n\t\t\t\t";
					echo $project."</font>";
					echo "<br>\n\t\t\t\t";
				}
				?></td></tr></table>

				<script type="text/javascript">
				// Javascript which needs to be run.
				var userProjectCount   = "<?php echo $userProjectCount; ?>";
				var systemProjectCount = "<?php echo $systemProjectCount; ?>";
				// Keep projects tab looking consitent during navigation.
				for (var key=1; key<(systemProjectCount+userProjectCount); key++) {
					if(document.getElementById("figure_"+key)){
						var show_button_element = document.getElementById("show_"+key);
						show_button_element.style.visibility='hidden';
					}
				}
				// build iframes for data loading and processing.
				<?php
				for ($key=0; $key<$userProjectCount; $key++) {
					$project  = $projectFolders[$key];
					$handle   = fopen($directory."users/".$user."/projects/".$project."/dataType.txt", "r");
					$dataType = fgets($handle);
					fclose($handle);

					// Generate entries for projects which are being processed.
					echo "    if (document.getElementById('frameContainer.p2_".$key."')) {\n\t\t";   // processing file.
					echo "        var show_button_element = document.getElementById('show_".$key."');\n\t\t";
					echo "            show_button_element.style.visibility='hidden';\n\t\t";
					echo "        var el_p_".$key." = document.getElementById('frameContainer.p2_".$key."');\n\t\t";
					echo "        el_p_".$key.".innerHTML='<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" style=\"height:20px\" src=\"project.working.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
					echo "        top.frames['p_".$key."'].user    = \"".$user."\";\n\t\t";
					echo "        top.frames['p_".$key."'].project = \"".$project."\";\n\t\t";
					echo "        top.frames['p_".$key."'].key     = \"p_".$key."\";\n\t";
					echo "    }\n\t\t";

					// Generate entries for projects which need data files loaded.
					echo "    if (document.getElementById('frameContainer.p3_".$key."')) {\n\t\t";   // starting... waiting for file upload
					echo "        var show_button_element = document.getElementById('show_".$key."');\n\t\t";
					echo "            show_button_element.style.visibility='hidden';\n\t\t";
					echo "        var el_p_".$key." = document.getElementById('frameContainer.p3_".$key."');\n\t\t";
					if ($dataType == '0') {
						echo "        // 0   = SnpCgh microarray.\n\t\t";
						echo "        el_p_".$key.".innerHTML='<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" style=\"height:38px\" src=\"uploader.1.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
						echo "        top.frames['p_".$key."'].display_string = new Array();\n\t\t";
						echo "        top.frames['p_".$key."'].display_string[0] = \"Add : SnpCgh array data...\";\n\t\t";
						echo "        top.frames['p_".$key."'].target_dir        = \"".$directory."users/".$user."/projects/".$project."/\";\n\t\t";
						echo "        top.frames['p_".$key."'].conclusion_script = \"".$url."php/project.SnpCgh.install.php\";\n\t\t";
						echo "        top.frames['p_".$key."'].user              = \"".$user."\";\n\t\t";
						echo "        top.frames['p_".$key."'].project           = \"".$project."\";\n\t\t";
						echo "        top.frames['p_".$key."'].key               = \"p_".$key."\";\n\t";

					} else if ($dataType == '1:0') {
						echo "        // 1:0 = WGseq : single-end [FASTQ/ZIP/GZ]\n\t\t";
						echo "        el_p_".$key.".innerHTML='<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" style=\"height:38px\" src=\"uploader.1.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
						echo "        top.frames['p_".$key."'].display_string = new Array();\n\t\t";
						echo "        top.frames['p_".$key."'].display_string[0] = \"Add : Single-end-read WGseq data (FASTQ/ZIP/GZ)...\";\n\t\t";
						echo "        top.frames['p_".$key."'].target_dir        = \"".$directory."users/".$user."/projects/".$project."/\";\n\t\t";
						echo "        top.frames['p_".$key."'].conclusion_script = \"".$url."php/project.single_WGseq.install_1.php\";\n\t\t";
						echo "        top.frames['p_".$key."'].user              = \"".$user."\";\n\t\t";
						echo "        top.frames['p_".$key."'].project           = \"".$project."\";\n\t\t";
						echo "        top.frames['p_".$key."'].key               = \"p_".$key."\";\n\t";
					} else if ($dataType == '1:1') {
						echo "        // 1:1 = WGseq : paired-end [FASTQ/ZIP/GZ]\n\t\t";
						echo "        el_p_".$key.".innerHTML='<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" style=\"height:76px\" src=\"uploader.2.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
						echo "        top.frames['p_".$key."'].display_string = new Array();\n\t\t";
						echo "        top.frames['p_".$key."'].display_string[0] = \"Add : Paired-end-read WGseq data (1/2; FASTQ/ZIP/GZ)...\";\n\t\t";
						echo "        top.frames['p_".$key."'].display_string[1] = \"Add : Paired-end-read WGseq data (2/2; FASTQ/ZIP/GZ)...\";\n\t\t";
						echo "        top.frames['p_".$key."'].target_dir        = \"".$directory."users/".$user."/projects/".$project."/\";\n\t\t";
						echo "        top.frames['p_".$key."'].conclusion_script = \"".$url."php/project.paired_WGseq.install_1.php\";\n\t\t";
						echo "        top.frames['p_".$key."'].user              = \"".$user."\";\n\t\t";
						echo "        top.frames['p_".$key."'].project           = \"".$project."\";\n\t\t";
						echo "        top.frames['p_".$key."'].key               = \"p_".$key."\";\n\t";
					} else if (($dataType == '1:2') || ($dataType == '1:3')) {
						echo "        // 1:2 or 1:3 = WGseq : [SAM/BAM/TXT]\n\t\t";
						echo "        el_p_".$key.".innerHTML='<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" style=\"height:38px\" src=\"uploader.1.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
						echo "        top.frames['p_".$key."'].display_string = new Array();\n\t\t";
						echo "        top.frames['p_".$key."'].display_string[0] = \"Add : WGseq data (SAM/BAM/TXT)...\";\n\t\t";
						echo "        top.frames['p_".$key."'].target_dir        = \"".$directory."users/".$user."/projects/".$project."/\";\n\t\t";
						echo "        top.frames['p_".$key."'].conclusion_script = \"".$url."php/project.single_WGseq.install_1.php\";\n\t\t";
						echo "        top.frames['p_".$key."'].user              = \"".$user."\";\n\t\t";
						echo "        top.frames['p_".$key."'].project           = \"".$project."\";\n\t\t";
						echo "        top.frames['p_".$key."'].key               = \"p_".$key."\";\n\t";

					} else if ($dataType == '2:0') {
						echo "        // 2:0 = ddRADseq : single-end [FASTQ/ZIP/GZ]\n\t\t";
						echo "        el_p_".$key.".innerHTML='<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" style=\"height:38px\" src=\"uploader.1.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
						echo "        top.frames['p_".$key."'].display_string = new Array();\n\t\t";
						echo "        top.frames['p_".$key."'].display_string[0] = \"Add : Single-end-read ddRADseq data (FASTQ/ZIP/GZ)...\";\n\t\t";
						echo "        top.frames['p_".$key."'].target_dir        = \"".$directory."users/".$user."/projects/".$project."/\";\n\t\t";
						echo "        top.frames['p_".$key."'].conclusion_script = \"".$url."php/project.single_ddRADseq.install_1.php\";\n\t\t";
						echo "        top.frames['p_".$key."'].user              = \"".$user."\";\n\t\t";
						echo "        top.frames['p_".$key."'].project           = \"".$project."\";\n\t\t";
						echo "        top.frames['p_".$key."'].key               = \"p_".$key."\";\n\t";
					} else if ($dataType == '2:1') {
						echo "        // 2:1 = ddRADseq : paired-end [FASTQ/ZIP/GZ]\n\t\t";
						echo "        el_p_".$key.".innerHTML='<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" style=\"height:76px\" src=\"uploader.2.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
						echo "        top.frames['p_".$key."'].display_string = new Array();\n\t\t";
						echo "        top.frames['p_".$key."'].display_string[0] = \"Add : Paired-end-read ddRADseq data (1/2; FASTQ/ZIP/GZ)...\";\n\t\t";
						echo "        top.frames['p_".$key."'].display_string[1] = \"Add : Paired-end-read ddRADseq data (2/2; FASTQ/ZIP/GZ)...\";\n\t\t";
						echo "        top.frames['p_".$key."'].target_dir        = \"".$directory."users/".$user."/projects/".$project."/\";\n\t\t";
						echo "        top.frames['p_".$key."'].conclusion_script = \"".$url."php/project.paired_ddRADseq.install_1.php\";\n\t\t";
						echo "        top.frames['p_".$key."'].user              = \"".$user."\";\n\t\t";
						echo "        top.frames['p_".$key."'].project           = \"".$project."\";\n\t\t";
						echo "        top.frames['p_".$key."'].key               = \"p_".$key."\";\n\t";
					} else if (($dataType == '2:2') || ($dataType == '2:3')) {
						echo "        // 2:2 or 2:3 = ddRADseq : [SAM/BAM/TXT]\n\t\t";
						echo "        el_p_".$key.".innerHTML='<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" style=\"height:38px\" src=\"uploader.1.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
						echo "        top.frames['p_".$key."'].display_string = new Array();\n\t\t";
						echo "        top.frames['p_".$key."'].display_string[0] = \"Add : ddRADseq data (SAM/BAM/TXT)...\";\n\t\t";
						echo "        top.frames['p_".$key."'].target_dir        = \"".$directory."users/".$user."/projects/".$project."/\";\n\t\t";
						echo "        top.frames['p_".$key."'].conclusion_script = \"".$url."php/project.single_ddRADseq.install_1.php\";\n\t\t";
						echo "        top.frames['p_".$key."'].user              = \"".$user."\";\n\t\t";
						echo "        top.frames['p_".$key."'].project           = \"".$project."\";\n\t\t";
						echo "        top.frames['p_".$key."'].key               = \"p_".$key."\";\n\t";

					} else if ($dataType == '3:0') {
                        echo "        // 3:0 = RNAseq : single-end [FASTQ/ZIP/GZ]\n\t\t";
                        echo "        el_p_".$key.".innerHTML='<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" style=\"height:38px\" src=\"uploader.1.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
                        echo "        top.frames['p_".$key."'].display_string = new Array();\n\t\t";
                        echo "        top.frames['p_".$key."'].display_string[0] = \"Add : Single-end-read RNAseq data (FASTQ/ZIP/GZ)...\";\n\t\t";
                        echo "        top.frames['p_".$key."'].target_dir        = \"".$directory."users/".$user."/projects/".$project."/\";\n\t\t";
                        echo "        top.frames['p_".$key."'].conclusion_script = \"".$url."php/project.single_RNAseq.install_1.php\";\n\t\t";
                        echo "        top.frames['p_".$key."'].user              = \"".$user."\";\n\t\t";
                        echo "        top.frames['p_".$key."'].project           = \"".$project."\";\n\t\t";
                        echo "        top.frames['p_".$key."'].key               = \"p_".$key."\";\n\t";
                    } else if ($dataType == '3:1') {
                        echo "        // 3:1 = RNAseq : paired-end [FASTQ/ZIP/GZ]\n\t\t";
                        echo "        el_p_".$key.".innerHTML='<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" style=\"height:76px\" src=\"uploader.2.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
                        echo "        top.frames['p_".$key."'].display_string = new Array();\n\t\t";
                        echo "        top.frames['p_".$key."'].display_string[0] = \"Add : Paired-end-read RNAseq data (1/2; FASTQ/ZIP/GZ)...\";\n\t\t";
                        echo "        top.frames['p_".$key."'].display_string[1] = \"Add : Paired-end-read RNAseq data (2/2; FASTQ/ZIP/GZ)...\";\n\t\t";
                        echo "        top.frames['p_".$key."'].target_dir        = \"".$directory."users/".$user."/projects/".$project."/\";\n\t\t";
                        echo "        top.frames['p_".$key."'].conclusion_script = \"".$url."php/project.paired_RNAseq.install_1.php\";\n\t\t";
                        echo "        top.frames['p_".$key."'].user              = \"".$user."\";\n\t\t";
                        echo "        top.frames['p_".$key."'].project           = \"".$project."\";\n\t\t";
                        echo "        top.frames['p_".$key."'].key               = \"p_".$key."\";\n\t";
                    } else if (($dataType == '3:2') || ($dataType == '3:3')) {
                        echo "        // 3:2 or 3:3 = RNAseq : [SAM/BAM/TXT]\n\t\t";
                        echo "        el_p_".$key.".innerHTML='<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" style=\"height:38px\" src=\"uploader.1.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
                        echo "        top.frames['p_".$key."'].display_string = new Array();\n\t\t";
                        echo "        top.frames['p_".$key."'].display_string[0] = \"Add : RNAseq data (SAM/BAM/TXT)...\";\n\t\t";
                        echo "        top.frames['p_".$key."'].target_dir        = \"".$directory."users/".$user."/projects/".$project."/\";\n\t\t";
                        echo "        top.frames['p_".$key."'].conclusion_script = \"".$url."php/project.single_RNAseq.install_1.php\";\n\t\t";
                        echo "        top.frames['p_".$key."'].user              = \"".$user."\";\n\t\t";
                        echo "        top.frames['p_".$key."'].project           = \"".$project."\";\n\t\t";
                        echo "        top.frames['p_".$key."'].key               = \"p_".$key."\";\n\t";

					} else if ($dataType == '4:0') {
						echo "        // 4:0 = IonExpressSeq : single-end [FASTQ/ZIP/GZ]\n\t\t";
						echo "        el_p_".$key.".innerHTML='<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" style=\"height:38px\" src=\"uploader.1.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
						echo "        top.frames['p_".$key."'].display_string = new Array();\n\t\t";
						echo "        top.frames['p_".$key."'].display_string[0] = \"Add : Single-end-read IonExpress data (FASTQ/ZIP/GZ)...\";\n\t\t";
						echo "        top.frames['p_".$key."'].target_dir        = \"".$directory."users/".$user."/projects/".$project."/\";\n\t\t";
						echo "        top.frames['p_".$key."'].conclusion_script = \"".$url."php/project.single_IonExpressSeq.install_1.php\";\n\t\t";
						echo "        top.frames['p_".$key."'].user              = \"".$user."\";\n\t\t";
						echo "        top.frames['p_".$key."'].project           = \"".$project."\";\n\t\t";
						echo "        top.frames['p_".$key."'].key               = \"p_".$key."\";\n\t";
					} else if ($dataType == '4:1') {
						echo "        // 4:1 = IonExpressSeq : paired-end [FASTQ/ZIP/GZ]\n\t\t";
						echo "        el_p_".$key.".innerHTML='<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" style=\"height:76px\" src=\"uploader.2.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
						echo "        top.frames['p_".$key."'].display_string = new Array();\n\t\t";
						echo "        top.frames['p_".$key."'].display_string[0] = \"Add : Paired-end-read IonExpress data (1/2; FASTQ/ZIP/GZ)...\";\n\t\t";
						echo "        top.frames['p_".$key."'].display_string[1] = \"Add : Paired-end-read IonExpress data (2/2; FASTQ/ZIP/GZ)...\";\n\t\t";
						echo "        top.frames['p_".$key."'].target_dir        = \"".$directory."users/".$user."/projects/".$project."/\";\n\t\t";
						echo "        top.frames['p_".$key."'].conclusion_script = \"".$url."php/project.paired_IonExpressSeq.install_1.php\";\n\t\t";
						echo "        top.frames['p_".$key."'].user              = \"".$user."\";\n\t\t";
						echo "        top.frames['p_".$key."'].project           = \"".$project."\";\n\t\t";
						echo "        top.frames['p_".$key."'].key               = \"p_".$key."\";\n\t";
					} else if (($dataType == '4:2') || ($dataType == '4:3')) {
						echo "        // 4:2 or 4:3 = IonExpressSeq : [SAM/BAM/TXT]\n\t\t";
						echo "        el_p_".$key.".innerHTML='<iframe id=\"p_".$key."\" name=\"p_".$key."\" class=\"upload\" style=\"height:38px\" src=\"uploader.1.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
						echo "        top.frames['p_".$key."'].display_string = new Array();\n\t\t";
						echo "        top.frames['p_".$key."'].display_string[0] = \"Add : IonExpress data (SAM/BAM/TXT)...\";\n\t\t";
						echo "        top.frames['p_".$key."'].target_dir        = \"".$directory."users/".$user."/projects/".$project."/\";\n\t\t";
						echo "        top.frames['p_".$key."'].conclusion_script = \"".$url."php/project.single_IonExpressSeq.install_1.php\";\n\t\t";
						echo "        top.frames['p_".$key."'].user              = \"".$user."\";\n\t\t";
						echo "        top.frames['p_".$key."'].project           = \"".$project."\";\n\t\t";
						echo "        top.frames['p_".$key."'].key               = \"p_".$key."\";\n\t";

					} else {
						// dataTypes not dealt with yet:
						//     4:0 = IonTorrent : single-end [FASTQ/ZIP/GZ]
						//     4:1 = IonTorrent : paired-end [FASTQ/ZIP/GZ]
						//     4:2 = IonTorrent : [SAM/BAM]
						//     other = unknown.
					}
					echo "} else {\n\t\t";   // complete.
					echo "}\n\t\t";
				}?>
				</script>
			</div>
<!---------------------------------------------------------------------------!>
			<div id="panel_genome" name="panel_genome">
				<table width="100%" cellpadding="0"><tr><td width="50%">
				<?php
				$userGenomeCount = 0;
				if (isset($_SESSION['logged_on'])) {
					$genomesDir    = $directory."users/".$user."/genomes/";
					$genomeFolders = array_diff(glob($genomesDir."*"), array('..', '.'));
					// Sort directories by date, newest first.
					array_multisort(array_map('filemtime', $genomeFolders), SORT_DESC, $genomeFolders);
					// Trim path from each folder string.
					foreach($genomeFolders as $key=>$folder) {   $genomeFolders[$key] = str_replace($genomesDir,"",$folder);   }
					$userGenomeCount = count($genomeFolders);
					echo "<button type='button' onclick=\"window.location.href='genome.create.php'\">Install New Genome</button><br>\n\t\t\t\t";
					echo "<b><font size='2'>User installed genomes:</font></b><br>\n\t\t\t\t";
					echo "<div class='tab' style='height:141px; overflow-y:scroll;'>\n\t\t\t\t";
					foreach($genomeFolders as $key=>$genome) {
						$genomeNameString = file_get_contents($directory."users/".$user."/genomes/".$genome."/name.txt");
						$genomeNameString = trim($genomeNameString);
						if (file_exists($directory."users/".$user."/genomes/".$genome."/complete.txt")) {
							echo "<span id='hapmap_label_".$key."' style='color:#00AA00;'>";
						} else if (file_exists($directory."users/".$user."/genomes/".$genome."/working.txt")) {
							echo "<span id='hapmap_label_".$key."' style='color:#BB9900;'>";
						} else {
							echo "<span id='hapmap_label_".$key."' style='color:#CC0000;'>";
						}
						echo "<font size='2'>".$genomeNameString."</font>";
						echo "<span id='g1_".$key."_delete'><button type='button' onclick=\"deleteGenomeConfirmation('".$genome."','g1_".$key."_delete')\"><font size='2'>Delete</font></button></span>\n\t\t\t\t";
						echo "<span id='g2_".$genome."_delete'></span><span id='g_".$genome."_type'></span><br></span>\n\t\t\t\t";
						if ((file_exists($directory."users/".$user."/genomes/".$genome."/complete.txt")) || (file_exists($directory."users/".$user."/genomes/".$genome."/working.txt"))) {
							echo "<div id='frameContainer.g1_".$key."'></div>\n\t\t\t\t";
						} else {
							echo "<div id='frameContainer.g2_".$key."'></div>\n\t\t\t\t";
						}
					}
				}?></div></td><td width="50%"><br><?php
				$genomesDir          = $directory."users/default/genomes/";
				$systemGenomeFolders = array_diff(glob($genomesDir."*"), array('..', '.'));
				// Sort directories by date, newest first.
				array_multisort(array_map('filemtime', $systemGenomeFolders), SORT_DESC, $systemGenomeFolders);
				// Trim path from each folder string.
				foreach($systemGenomeFolders as $key=>$folder) {   $systemGenomeFolders[$key] = str_replace($genomesDir,"",$folder);   }
				$systemGenomeCount = count($systemGenomeFolders);
				echo "<b><font size='2'>System genomes:</font></b><br>";
				echo "<div class='tab' style='height:141px; overflow-y:scroll;'>";
				foreach ($systemGenomeFolders as $key=>$genome) {
					$genomeNameString = file_get_contents($directory."users/default/genomes/".$genome."/name.txt");
					$genomeNameString = trim($genomeNameString);
					echo "<font size='2'>".$genomeNameString."</font><br>";
				}
				?></td></tr></table>

				<script type="text/javascript">
				// Javascript which needs to be run.
				var userGenomeCount   = "<?php echo $userGenomeCount; ?>";
				var systemGenomeCount = "<?php echo $systemGenomeCount; ?>";
				<?php
				for ($key=0; $key<$userGenomeCount; $key++) {
					$genome = $genomeFolders[$key];
					echo "if (document.getElementById('frameContainer.g2_".$key."')) {\n\t\t";   // waiting for file upload.
					echo "    var el_g_".$key." = document.getElementById('frameContainer.g2_".$key."');\n\t\t";
					echo "    el_g_".$key.".innerHTML='<iframe id=\"g_".$key."\" name=\"g_".$key."\" class=\"upload\" style=\"height:38px\" src=\"uploader.1.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
					echo "    top.frames['g_".$key."'].display_string    = new Array();\n\t\t";
					echo "    top.frames['g_".$key."'].display_string[0] = \"Add : Genome reference FASTA file...\";\n\t\t";
					echo "    top.frames['g_".$key."'].target_dir        = \"".$directory."users/".$user."/genomes/".$genome."/\";\n\t\t";
					echo "    top.frames['g_".$key."'].conclusion_script = \"".$url."php/genome.install_1.php\";\n\t\t";
					echo "} else {\n\t\t";   // working or complete.
					echo "    var el_g_".$key." = document.getElementById(\"frameContainer.g1_".$key."\");\n\t\t";
					echo "    el_g_".$key.".innerHTML='<iframe id=\"g_".$key."\" name=\"g_".$key."\" class=\"upload\" src=\"genome.working.php\" marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t";
					echo "}\n\t\t";
					echo "top.frames['g_".$key."'].user              = \"".$user."\";\n\t\t";
					echo "top.frames['g_".$key."'].genome            = \"".$genome."\";\n\t\t";
					echo "top.frames['g_".$key."'].key               = \"g_".$key."\";\n\t";
				}
				?>
				</script>
			</div>
<!---------------------------------------------------------------------------!>
			<div id="panel_hapmap" name="panel_hapmap">
				<script type="text/javascript">
				function showColors(colorName,targetToChange,contentString) {
					if(document.getElementById(targetToChange)){
						var select       = document.getElementById(targetToChange);
						select.innerHTML = '';
						fillerString     = '<div style=\"width:35px; height:20px; overflow:hidden; display:inline-block; background-color:';
						if        (colorName == "deep pink")       { fillerString = fillerString + 'rgb(255,0  ,127)';
						} else if (colorName == "magenta")         { fillerString = fillerString + 'rgb(255,0  ,255)';
						} else if (colorName == "electric indigo") { fillerString = fillerString + 'rgb(127,0  ,255)';
						} else if (colorName == "blue")            { fillerString = fillerString + 'rgb(0  ,0  ,255)';
						} else if (colorName == "dodger blue")     { fillerString = fillerString + 'rgb(0  ,127,255)';
						} else if (colorName == "cyan")            { fillerString = fillerString + 'rgb(0  ,255,255)';
						} else if (colorName == "spring green")    { fillerString = fillerString + 'rgb(0  ,255,127)';
						} else if (colorName == "green")           { fillerString = fillerString + 'rgb(0  ,255,0  )';
						} else if (colorName == "chartreuse")      { fillerString = fillerString + 'rgb(127,255,0  )';
						} else if (colorName == "yellow")          { fillerString = fillerString + 'rgb(255,255,0  )';
						} else if (colorName == "dark orange")     { fillerString = fillerString + 'rgb(255,127,0  )';
						} else if (colorName == "red")             { fillerString = fillerString + 'rgb(255,0  ,0  )';
						} else { fillerString = fillerString + 'rgb(127,127,127)';
						}
						fillerString     = fillerString + '\"><center><b>' + contentString + '</b></center></div>';
						select.innerHTML = fillerString;
					}
				}
				</script>
				<table width="100%" cellpadding="0"><tr valign="top"><td width="50%">
				<?php
				if (isset($_SESSION['user']) != 0) {
					$hapmapsDir    = $directory."users/".$user."/hapmaps/";
					$hapmapFolders = array_diff(glob($hapmapsDir."*"), array('..', '.'));
					// Sort directories by date, newest first.
					array_multisort(array_map('filemtime', $hapmapFolders), SORT_DESC, $hapmapFolders);
					// Trim path from each folder string.
					foreach($hapmapFolders as $key=>$folder) {   $hapmapFolders[$key] = str_replace($hapmapsDir,"",$folder);   }
					echo "<div class='hapmap'><b><font size='2'>User generated hapmaps:</font></b><br>\n\t\t\t\t";
					echo "<button class='tab' type='button' onclick=\"window.location.href='hapmap.create.php'\">Generate New Hapmap</button>\n\t\t\t\t";
					foreach($hapmapFolders as $key=>$hapmap) {
						echo "<div class='tab'><table><tr><td>\n\t\t\t\t";
						if (file_exists($directory."users/".$user."/hapmaps/".$hapmap."/complete.txt")) {
							echo "<span id='hapmap_label_".$key."' style='color:#00AA00;'>";
						} else if (file_exists($directory."users/".$user."/hapmaps/".$hapmap."/working.txt")) {
							echo "<span id='hapmap_label_".$key."' style='color:#BB9900;'>";
						} else {
							echo "<span id='hapmap_label_".$key."' style='color:#CC0000;'>";
						}
						echo "<font size='2'>".$hapmap."</font></span></td><td>\n\t\t\t\t";
						echo "<span id='h1_".$key."_delete'><button type='button' onclick=\"deleteHapmapConfirmation('".$hapmap."','h1_".$key."_delete')\">Delete</button></span>\n\t\t\t\t";
						echo "<span id='h2_".$hapmap."_delete'></span></td>\n\t\t\t\t";
						echo "<td><div id='userHapmapA_".$key."'  >filler</div></td>\n\t\t\t\t";
						echo "<td><div id='userHapmapHET_".$key."'>filler</div></td>\n\t\t\t\t";
						echo "<td><div id='userHapmapB_".$key."'  >filler</div></td>\n\t\t\t\t";
						echo "<td><div id='userHapmapHOM_".$key."'>filler</div></td></tr></table>\n\t\t\t\t";
						echo "<div id='frameContainer.h_".$key."'></div></div>\n\t\t\t\t";

						echo "<script type='text/javascript'>\n\t\t\t\t";
						echo "var el_h_".$key." = document.getElementById('frameContainer.h_".$key."');\n\t\t\t\t";
						echo "function hapmap_UI_refresh_1_".$key."() {\n\t\t\t\t";
						echo "    var frameString3 = '<iframe id=\"h_".$key."\" name=\"h_".$key."\" class=\"upload\" style=\"height:38px\" src=\"hapmap.addTo.php\"';\n\t\t\t\t";
						echo "    var frameString4 = ' marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t\t\t";
						echo "    var frameString5 = '<iframe id=\"h2_".$key."\" name=\"h2_".$key."\" class=\"upload\" style=\"height:38px\" src=\"hapmap.finalize.php\"';\n\t\t\t\t";
						echo "    var frameString6 = '  marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t\t\t";
						echo "    el_h_".$key.".innerHTML = frameString3.concat(frameString4,frameString5,frameString6);\n\t\t\t\t";
						echo "    top.frames['h_".$key."'].hapmap            = '".$hapmap."';\n\t\t\t\t";
						echo "    top.frames['h_".$key."'].user              = '".$user."';\n\t\t\t\t";
						echo "    top.frames['h_".$key."'].key               = 'h_".$key."';\n\t\t\t\t";
						echo "    top.frames['h2_".$key."'].hapmap           = '".$hapmap."';\n\t\t\t\t";
						echo "    top.frames['h2_".$key."'].user             = '".$user."';\n\t\t\t\t";
						echo "    top.frames['h2_".$key."'].key              = 'h2_".$key."';\n\t\t\t\t";
						echo "}\n\t\t\t\t";
						echo "function hapmap_UI_refresh_2_".$key."() {\n\t\t\t\t";
						echo "    var frameString1 = '<iframe id=\"h_".$key."\" name=\"h_".$key."\" class=\"upload\" src=\"hapmap.working.php\" ';\n\t\t\t\t";
						echo "    var frameString2 = 'marginwidth=\"0\" marginheight=\"0\" vspace=\"0\" hspace=\"0\"></iframe>';\n\t\t\t\t";
						echo "    el_h_".$key.".innerHTML = frameString1.concat(frameString2);\n\t\t\t\t";
						echo "    top.frames['h_".$key."'].hapmap            = '".$hapmap."';\n\t\t\t\t";
						echo "    top.frames['h_".$key."'].user              = '".$user."';\n\t\t\t\t";
						echo "    top.frames['h_".$key."'].key               = 'h_".$key."';\n\t\t\t\t";
						echo "}\n\t\t\t\t";

						if (file_exists($directory."/users/".$user."/hapmaps/".$hapmap."/complete.txt")) {
							// No user interface required.
						} elseif (file_exists($directory."/users/".$user."/hapmaps/".$hapmap."/working.txt")) {
							// The hapmap is being processed.
							// Load 'hapmap.working.php' into an internal frame 'h_$key'.
							echo "hapmap_UI_refresh_2_".$key."();\n\t\t\t\t";
						} else {
							echo "hapmap_UI_refresh_1_".$key."();\n\t\t\t\t";
						}
						echo "</script>\n\t\t\t\t";
					}
				}
				if (file_exists($directory."/users/".$user."/hapmaps/".$hapmap."/complete.txt")) {
					// No user interface required.
				} elseif (file_exists($directory."/users/".$user."/hapmaps/".$hapmap."/working.txt")) {
					// The hapmap is being processed.
				} else {
					// hapmap addto function needed.
				}?></div></td><td width="50%"><br>
				<?php
				$hapmapsDir          = $directory."users/default/hapmaps/";
				$systemHapmapFolders = array_diff(glob($hapmapsDir."*"), array('..', '.'));
				// Sort directories by date, newest first.
				array_multisort(array_map('filemtime', $systemHapmapFolders), SORT_DESC, $systemHapmapFolders);
				// Trim path from each folder string.
				foreach($systemHapmapFolders as $key=>$folder) {   $systemHapmapFolders[$key] = str_replace($hapmapsDir,"",$folder);   }
				echo "<div class='hapmap'><b><font size='2'>System hapmaps:</font></b>\n\t\t\t\t";
				foreach ($systemHapmapFolders as $key=>$hapmap) {
					echo "<div class='tab'><table><tr><td><font size='2'>".$hapmap."</font></td>\n\t\t\t\t";
					echo "<td><div id='defaultHapmapA_".$key."'  >filler</div></td>\n\t\t\t\t";
					echo "<td><div id='defaultHapmapHET_".$key."'>filler</div></td>\n\t\t\t\t";
					echo "<td><div id='defaultHapmapB_".$key."'  >filler</div></td>\n\t\t\t\t";
					echo "<td><div id='defaultHapmapHOM_".$key."'>filler</div></td></tr></table>\n\t\t\t\t";
					echo "</div>\n\t\t\t\t";
				}
				?></td></tr></table>

				<script type="text/javascript">
				// javascript which needs to run for the hapmap panel.
				<?php
				if (isset($_SESSION['user']) != 0) {
					$hapmapsDir    = $directory."users/".$user."/hapmaps/";
					$hapmapFolders = array_diff(glob($hapmapsDir."*"), array('..', '.'));
					// Sort directories by date, newest first.
					array_multisort(array_map('filemtime', $hapmapFolders), SORT_DESC, $hapmapFolders);
					// Trim path from each folder string.
					foreach($hapmapFolders as $key=>$folder) {   $hapmapFolders[$key] = str_replace($hapmapsDir,"",$folder);   }
					foreach($hapmapFolders as $key=>$hapmap) {
						$handle       = fopen($directory."users/".$user."/hapmaps/".$hapmap."/colors.txt",'r');
						$colorString1 = trim(fgets($handle));
						$colorString2 = trim(fgets($handle));
						fclose($handle);
						echo "\tshowColors('".$colorString1."','userHapmapA_".$key."','a');\n";
						echo "\tshowColors('".$colorString2."','userHapmapB_".$key."','b');\n";
						echo "\tshowColors('red','userHapmapHOM_".$key."','hom');\n";
						echo "\tshowColors('grey','userHapmapHET_".$key."','ab');\n";
					}
				}
				$hapmapsDir          = $directory."users/default/hapmaps/";
				$systemHapmapFolders = array_diff(glob($hapmapsDir."*"), array('..', '.'));
				// Sort directories by date, newest first.
				array_multisort(array_map('filemtime', $systemHapmapFolders), SORT_DESC, $systemHapmapFolders);
				// Trim path from each folder string.
				foreach($systemHapmapFolders as $key=>$folder) {   $systemHapmapFolders[$key] = str_replace($hapmapsDir,"",$folder);   }
				foreach ($systemHapmapFolders as $key=>$hapmap) {
					$handle       = fopen($directory."users/default/hapmaps/".$hapmap."/colors.txt",'r');
					$colorString1 = trim(fgets($handle));
					$colorString2 = trim(fgets($handle));
					fclose($handle);
					echo "\tshowColors('".$colorString1."','defaultHapmapA_".$key."','a');\n";
					echo "\tshowColors('".$colorString2."','defaultHapmapB_".$key."','b');\n";
					echo "\tshowColors('red','defaultHapmapHOM_".$key."','hom');\n";
					echo "\tshowColors('grey','defaultHapmapHET_".$key."','ab');\n";
				}?>
				</script>
			</div>
<!---------------------------------------------------------------------------!>
			<div id="panel_system" name="panel_system">
				<b>Interact with an admin:</b><br>
				<div class="tab">
				<?php
				if (isset($_SESSION['user']) != 0) {
					echo "Report a bug! Request a feature!<br>";
					echo "Or just as the admin a question about how to use a certain feature.<br>";
					echo "<button type='button' onclick=\"window.location.href='bug.php'\">Click here...</button>";
				} else {
					echo "Log in using the 'User' tab to gain access to the admin.";
				}?>
				</div><br>
				<b>System status:</b><br><div id="frameContainer.status"></div>
				<script type="text/javascript">
				var el_status = document.getElementById("frameContainer.status");
				<?php
					// Load last line from "status_log.txt" file.
					$statusLogLines = explode("\n", trim(file_get_contents($directory."status_log.txt")));
					$statusLog      = str_replace("\n","<br>",trim(file_get_contents($directory."status_log.txt")));
					if (count($statusLogLines) > 1) {
						$statusLogEntry = $statusLog;
					} else if (count($statusLogLines) == 1) {
						if ($statusLog[count($statusLog)-1] == "") {
							$statusLogEntry = "No system status notes.";
						} else {
							$statusLogEntry = $statusLog;
						}
					} else {
						$statusLogEntry = "No system status notes.";
					}
				?>el_status.innerHTML = '<div class="tab"><font size="2"><?php echo $statusLogEntry; ?></font></div>';
				</script>
			</div>
<!---------------------------------------------------------------------------!>
			<div id="panel_about" name="panel_about">
			<div class='tab' style='height:185px; overflow-y:scroll;'>
				<font size="2">
				<p>The Yeast Mapping Analysis Pipeline (YMAP) is intended to simplify the processing of whole genome array and sequence
				datasets from yeast (or other small genome organisms) into simple to interpret figures illustrating large-scale genomic
				alterations.</p>
				<p>The current version (v1.0) of the pipeline is designed to process the following types of data:
				<ol>
					<li>SnpCgh Microarray.</li>
					<li>Whole genome next-generation-seq (WGseq).</li>
					<li>Double digest restriction site associated sequencing (ddRADseq).</li>
				</ol>
				</p>
				<hr width="50%">
				<p>The input file format options for each data type are as follows :</p>
				<p><b>1. SnpCgh microarray</b>
					<ul>
						<li>Tab-delimited text with the following format :
						<table width="600px" class="tab"><tr><td>
							1.&nbsp;Data starts in the 47th row.<br>
							2.&nbsp;1st column holds probe ID/name.<br>
							3.&nbsp;4th column holds probe channel 1 data.<br>
							4.&nbsp;5th column holds probe channel 2 data.<br>
							5.&nbsp;6th column holds probe channel ratio data.<br>
							6.&nbsp;7th column holds probe channel log2ratio data.
						</td></tr></table>
						</li>
					</ul>
				</p>
				<p><b>2. WGseq</b> and <b>3. ddRADseq</b>
					<ul>
						<li>Single-end reads as one file in raw FASTQ (*.fastq; *.fq) or compressed (*.zip; *.gz) format.</li>
						<li>Paired-end reads as two files in raw FASTQ (*.fastq; *.fq) or compressed (*.zip; *.gz) format.</li>
						<li>Sequence Alignment/Map file in raw (*.sam) or compressed (*.bam) format.</li>
						<li>For examination of custom-filtered data, tab-delimited text with the following format can be used :
						<table width="600px" class="tab"><tr><td>
							1.&nbsp;Chromosome name (must match reference).<br>
							2.&nbsp;Chromosome bp coordinate.<br>
							3.&nbsp;1' base-call.<br>
							4.&nbsp;Count of 1' base-call.<br>
							5.&nbsp;<i>2' base-call.</i>
						</td><td>
							6.&nbsp;<i>Count of 2' base-call.</i><br>
							7.&nbsp;<i>3' base-call.</i><br>
							8.&nbsp;<i>Count of 3' base-call.</i><br>
							9.&nbsp;<i>4' base-call.</i><br>
							10.&nbsp;<i>Count of 4' base-call.</i>
						</td></tr></table>
						</li>
					</ul>
				</p>
				<hr width="50%">
				<p>Projects are listed with those needing data file upload first, followed by those which are in-process, and then those that are complete. As projects transition from
				one category to the next, the color of their line will be updated. After a page refresh, position in the list will be updated.</p>
				<p>Genomes can be installed by the user. The tool has been designed for the analysis of relatively small (~2Gbase) fungal genomes, but others can be used.</p>
				<p>Haplotype maps can be constructed starting from one heterozygous parent or two homozygous parent datasets. When starting from a single heterozygous parent, additional
				child datasets with large-scale homozygoses are needed to construct the final map.</p>
			</font>
			</div>
			</div>
<!---------------------------------------------------------------------------!>
			<div id="panel_examples" name="panel_examples">
			<div class='tab' style='height:185px; overflow-y:scroll;'>
			<font size="2">
				<br>
				Datasets collected in association with the paper are linked below in the format of "<b>file:[file type]</b>".<br>
				You can read about file format requirements for Ymap by selecting the "About" tab above.<br><br>

				Some datasets were retrieved from external databases and are linked below in the format of "<b>url:[database name, accession]</b>".<br>
				To download datasets from the NCBI-SRA database in FASTQ format, enter the experiment number for the sample into <a href="http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=search_seq_name">[this page]</a>, then select the accession number of interest.<br><br>

				<b>Figure 2 : WGseq, chr-end bias.</b><br>
					<ul>
					<li>YQ2, url:<a href="http://www.ebi.ac.uk/biosamples/sample/SAMEA1879786">[EMBL-EBI BiosSamples, SAMEA1879786]</a>
						  or url:<a href="http://www.ncbi.nlm.nih.gov/sra/ERX238051">[NCBI-SRA, ERX238051]</a></li>
					</ul>
				<b>Figure 3 : WGseq, GC-content bias.</b><br>
					<ul>
					<li>FH6, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig3.FH6_forward.fastq.zip">[ZIP compressed FASTQ]</a>
						 and file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig3.FH6_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					</ul>
				<b>Figure 4 : ddRADseq.</b><br>
					<ul>
					<li>SC5314, file:<a   href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig4.SC5314.txt">[TXT]</a></li>
					<li>YJB11461, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig4.YJB11461.txt">[TXT]</a></li>
					<li>YJB9442, file:<a  href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig4.YJB9442.txt">[TXT]</a></li>
					</ul>
				<b>Figure 5 : WGseq</b><br>
					<ul>
					<li>SC5314, url:<a    href="http://www.ncbi.nlm.nih.gov/sra/?term=SRR868699">[NCBI-SRA, SRR868699]</a></li>
					<li>FH5, file:<a      href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig5.FH5_forward.fastq.zip">[ZIP compressed FASTQ]</a>
						 and file:<a      href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig5.FH5_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					<li>YJB12746, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig5.YJB12746_forward.fastq.zip">[ZIP compressed FASTQ]</a>
							  and file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig5.YJB12746_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					</ul>
				<b>Figure 6 : ddRADseq</b><br>
					<ul>
					<li>SC5314, file:<a                 href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig6.SC5314.txt">[TXT]</a></li>
					<li>YJB12712 derivative #1, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig6.YJB12712_derivative_1.txt">[TXT]</a></li>
					<li>YJB12712 derivative #2, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig6.YJB12712_derivative_2.txt">[TXT]</a></li>
					<li>YJB12712 derivative #9, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig6.YJB12712_derivative_9.txt">[TXT]</a></li>
					</ul>
				<b>Figure 8 : comparisons between array and sequence datasets</b><br>
					<ul>
					<li>YJB10490 array, file:<a    href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig8.YJB10490_array.xls">[XLS]</a></li>
					<li>YJB10490 WGseq, file:<a    href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig8.YJB10490_forward.fastq.zip">[ZIP compressed FASTQ]</a>
									and file:<a    href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig8.YJB10490_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					<li>YJB12229 array, file:<a    href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig8.YJB12229_array.xls">[XLS]</a></li>
					<li>YJB12229 ddRADseq, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig8.YJB12229.txt">[TXT]</a></li>
					<li>Ss2 array, file:<a         href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig8.Ss2_array.xls">[XLS]</a></li>
					<li>YJB12353 WGseq, file:<a    href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig8.YJB12353_forward.fastq.zip">[ZIP compressed FASTQ</a>
									and file:<a    href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig8.YJB12353_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					</ul>
				<b>Figure 9 : LOH in clinical isolates</b><br>
					<ul>
					<li>SC5314, url:<a     href="http://www.ebi.ac.uk/biosamples/sample/SAMN02141741">[EMBL-EBI BioSamples, SAMN02141741]</a></li>
					<li>SC5314, file:<a    href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig9.SC5314_forward.fastq.zip">[ZIP compressed FASTQ]</a>
							and file:<a    href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig9.SC5314_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					<li>SC5314, url:<a     href="http://www.ncbi.nlm.nih.gov/sra/?term=SRR868699">[NCBI-SRA, SRR868699]</a></li>
					<li>FH1, file:<a       href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig9.FH1_forward.fastq.zip">[ZIP compressed FASTQ]</a>
						 and file:<a       href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig9.FH1_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					<li>ATCC200955, url:<a href="http://www.ncbi.nlm.nih.gov/sra/?term=SAMN02140345">[NCBI-SRA, SAMN02140345]</a></li>
					<li>ATCC10231, url:<a  href="http://www.ncbi.nlm.nih.gov/sra/?term=SAMN02140347">[NCBI-SRA, SAMN02140347]</a></li>
					<li>YL1, url:<a        href="http://www.ebi.ac.uk/biosamples/sample/SAMEA1879767">[EMBL-EBI BioSamples, SAMEA1879767]</a></li>
					<li>YQ2, url:<a        href="http://www.ebi.ac.uk/biosamples/sample/SAMEA1879786">[EMBL-EBI BiosSamples, SAMEA1879786]</a>
						  or url:<a        href="http://www.ncbi.nlm.nih.gov/sra/ERX238051">[NCBI-SRA, ERX238051]</a></li>
					</ul>
				<b>Figure 10 : FH series</b>
					<ul>
					<li>FH1, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH1_forward.fastq.zip">[ZIP compressed FASTQ]</a>
						 and file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH1_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					<li>FH2, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH2_forward.fastq.zip">[ZIP compressed FASTQ]</a>
						 and file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH2_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					<li>FH3, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH3_forward.fastq.zip">[ZIP compressed FASTQ]</a>
						 and file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH3_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					<li>FH4, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH4_forward.fastq.zip">[ZIP compressed FASTQ]</a>
						 and file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH4_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					<li>FH5, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH5_forward.fastq.zip">[ZIP compressed FASTQ]</a>
						 and file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH5_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					<li>FH6, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH6_forward.fastq.zip">[ZIP compressed FASTQ]</a>
						 and file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH6_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					<li>FH7, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH7_forward.fastq.zip">[ZIP compressed FASTQ]</a>
						 and file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH7_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					<li>FH8, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH8_forward.fastq.zip">[ZIP compressed FASTQ]</a>
						 and file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH8_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					<li>FH9, file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH9_forward.fastq.zip">[ZIP compressed FASTQ]</a>
						 and file:<a href="http://lovelace.cs.umn.edu/Ymap/example_datasets/Fig10.FH9_reverse.fastq.zip">[ZIP compressed FASTQ]</a></li>
					</ul>
			</font>
			</div>
			</div>
<!---------------------------------------------------------------------------!>
			</td>
		</tr></table>
<!---------------------------------------------------------------------------!>
<!--- END : Construct tabbed menu and content pane.  ------------------------!>
<!---------------------------------------------------------------------------!>

	</td></tr></table>

	<script type="text/javascript">
	<!--------------- Tab management functions --------------!>
	deselect_tab = function(name) {
		var current_tab = document.getElementById("tab_"+name);
		    current_tab.style.border="1px solid #000000";
		    current_tab.style.backgroundColor="#DDDDDD";
		if (name != "user") {
			current_tab.style.borderLeft="none";
		}
		document.getElementById("panel_"+name).style.display='none';
	}
	select_tab = function(name) {
		var current_tab = document.getElementById("tab_"+name);
		    current_tab.style.border="1px solid #000000";
		    current_tab.style.borderBottom="none";
		    current_tab.style.backgroundColor="#FFFFFF";
		if (name != "user") {
			current_tab.style.borderLeft="none";
		}
		document.getElementById("panel_"+name).style.display='inline';
	}
	blank_and_content_tab = function() {
		var current_tab = document.getElementById("tab_blank");
		    current_tab.style.borderBottom="1px solid #000000";
		var current_tab = document.getElementById("tab_content");
		    current_tab.style.border="1px solid #000000";
		    current_tab.style.borderTop="none";
	}
	tabWindow_user = function() {
		select_tab("user");
		deselect_tab("project");
		deselect_tab("genome");
		deselect_tab("hapmap");
		deselect_tab("system");
		deselect_tab("about");
		deselect_tab("examples");
		blank_and_content_tab();
		var tabInUse = "user";
		localStorage.setItem("tabInUse", tabInUse);
	}
	tabWindow_project = function() {
		deselect_tab("user");
		select_tab("project");
		deselect_tab("genome");
		deselect_tab("hapmap");
		deselect_tab("system");
		deselect_tab("about");
		deselect_tab("examples");
		blank_and_content_tab();
		var tabInUse = "project";
		localStorage.setItem("tabInUse", tabInUse);
	}
	tabWindow_genome = function() {
		deselect_tab("user");
		deselect_tab("project");
		select_tab("genome");
		deselect_tab("hapmap");
		deselect_tab("system");
		deselect_tab("about");
		deselect_tab("examples");
		blank_and_content_tab();
		var tabInUse = "genome";
		localStorage.setItem("tabInUse", tabInUse);
	}
	tabWindow_hapmap = function() {
		deselect_tab("user");
		deselect_tab("project");
		deselect_tab("genome");
		select_tab("hapmap");
		deselect_tab("system");
		deselect_tab("about");
		deselect_tab("examples");
		blank_and_content_tab();
		var tabInUse = "hapmap";
		localStorage.setItem("tabInUse", tabInUse);
	}
	tabWindow_system = function() {
		deselect_tab("user");
		deselect_tab("project");
		deselect_tab("genome");
		deselect_tab("hapmap");
		select_tab("system");
		deselect_tab("about");
		deselect_tab("examples");
		blank_and_content_tab();
		var tabInUse = "system";
		localStorage.setItem("tabInUse", tabInUse);
	}
	tabWindow_about = function() {
		deselect_tab("user");
		deselect_tab("project");
		deselect_tab("genome");
		deselect_tab("hapmap");
		deselect_tab("system");
		select_tab("about");
		deselect_tab("examples");
		blank_and_content_tab();
		var tabInUse = "about";
		localStorage.setItem("tabInUse", tabInUse);
	}
	tabWindow_examples = function() {
		deselect_tab("user");
		deselect_tab("project");
		deselect_tab("genome");
		deselect_tab("hapmap");
		deselect_tab("system");
		deselect_tab("about");
		select_tab("examples");
		blank_and_content_tab();
		var tabInUse = "examples";
		localStorage.setItem("tabInUse", tabInUse);
	}

	<!---------------- Project display section --------------!>
	isFile = function(str) {
		var O= AJ();
		if(!O) return false;
		try {
			O.open("HEAD", str, false);
			O.send(null);
			return (O.status==200) ? true : false;
		}
		catch(er) { return false; }
	}
	AJ = function() {
		var obj;
		if (window.XMLHttpRequest){
			obj= new XMLHttpRequest();
		} else if (window.ActiveXObject) {
			try{       obj= new ActiveXObject('MSXML2.XMLHTTP.3.0'); }
			catch(er){ obj=false; }
		}
		return  obj;
	}
	loadImage = function(key,imageUrl,imageScale) {
		document.getElementById('fig_'+key).innerHTML = "<img src='"+imageUrl+"' width='"+imageScale+"%'>";
	}
	loadExternal = function(imageUrl) {
		window.open(imageUrl);
	}
	openProject = function(user,project,key,color1,color2,parent) {
		var url_base                         = "<?php echo $url; ?>";
		var fig_linear_CNV_SNP               = url_base+"users/"+user+"/projects/"+project+"/fig.CNV-SNP-map.2.";
		var fig_standard_CNV_SNP             = url_base+"users/"+user+"/projects/"+project+"/fig.CNV-SNP-map.1.";
		var fig_linear_CNV                   = url_base+"users/"+user+"/projects/"+project+"/fig.CNV-map.2.";
		var fig_standard_CNV                 = url_base+"users/"+user+"/projects/"+project+"/fig.CNV-map.1.";
		var fig_linear_manual                = url_base+"users/"+user+"/projects/"+project+"/fig.CNV-manualLOH-map.2.";
		var fig_standard_manual              = url_base+"users/"+user+"/projects/"+project+"/fig.CNV-manualLOH-map.1.";
		var CGD_annotations_SNP              = url_base+"users/"+user+"/projects/"+project+"/CGD_annotations."+project+".txt";

		var CNV_bias_WGseq_1                 = url_base+"users/"+user+"/projects/"+project+"/fig.examine_bias.png";
		var CNV_bias_WGseq_2                 = url_base+"users/"+user+"/projects/"+project+"/fig.GCratio_vs_CGH.png";

		var CNV_bias_ddRADseq_1              = url_base+"users/"+user+"/projects/"+project+"/fig.examine_bias.1.png";
		var CNV_bias_ddRADseq_2              = url_base+"users/"+user+"/projects/"+project+"/fig.examine_bias.2.png";
		var CNV_bias_ddRADseq_3              = url_base+"users/"+user+"/projects/"+project+"/fig.examine_bias.3.png";
		var CNV_bias_ddRADseq_4              = url_base+"users/"+user+"/projects/"+project+"/fig.examine_bias.T.png";

		var fig_a_linear_SNPratio_ddRADseq   = url_base+"users/"+user+"/projects/"+project+"/fig.allelic_fraction_histogram.png";
		var fig_b_linear_SNPratio_ddRADseq   = url_base+"users/"+user+"/projects/"+project+"/fig.allelic_ratio-map.b2.png";
		var fig_b_standard_SNPratio_ddRADseq = url_base+"users/"+user+"/projects/"+project+"/fig.allelic_ratio-map.b1.png";
		var fig_c_linear_SNPratio_ddRADseq   = url_base+"users/"+user+"/projects/"+project+"/fig.allelic_ratio-map.c2.png";
		var fig_c_standard_SNPratio_ddRADseq = url_base+"users/"+user+"/projects/"+project+"/fig.allelic_ratio-map.c1.png";
		var fig_d_linear_SNPratio_ddRADseq   = url_base+"users/"+user+"/projects/"+project+"/fig.allelic_ratio-map.d2.png";
		var fig_d_standard_SNPratio_ddRADseq = url_base+"users/"+user+"/projects/"+project+"/fig.allelic_ratio-map.d1.png";

		var visible_list         = document.getElementById("visible_list");
		if (isFile(fig_linear_CNV_SNP+"png")) {
			var mainFigure1 = fig_linear_CNV_SNP;
			var mainFigure2 = fig_standard_CNV_SNP;
		} else {
			var mainFigure1 = fig_linear_CNV;
			var mainFigure2 = fig_standard_CNV;
		}
		var string1 = "<div id='figure_"+key+"'><table border='0' align='center' width='100%'><tr><td width='35%' align='left'>";

		string1 = string1 + "<table><tr><td valign='top'>";
		string1 = string1 + project+" ";
		string1 = string1 + "</td><td valign='bottom'>";
		string1 = string1 + "<div id='userProjectA_"+key+"'   style='display:inline'></div>";
		string1 = string1 + "<div id='userProjectHET_"+key+"' style='display:inline'></div>";
		string1 = string1 + "<div id='userProjectB_"+key+"'   style='display:inline'></div>";
		string1 = string1 + "<div id='userProjectHOM_"+key+"' style='display:inline'></div>";
		string1 = string1 + "</td></tr></table>";

		string1 = string1 + "</td><td width='60%' align='left'>";
		string1 = string1 + "<font size='-1'>Linear</font> ";
		string1 = string1 + "<img src='images/icon_png_15b.png' alt-text='[PNG] button' align='center' onclick='loadImage(\""+key+"\",\""+mainFigure1+"png\",\"100\")'> ";
		string1 = string1 + "<img src='images/icon_eps_15b.png' alt-text='[EPS] button' align='center' onclick='loadExternal(\""+mainFigure1+"eps\")'>";
		string1 = string1 + " : <font size='-1'>Standard</font> ";
		string1 = string1 + "<img src='images/icon_png_15b.png' alt-text='[PNG] button' align='center' onclick='loadImage(\""+key+"\",\""+mainFigure2+"png\",\"50\")'> ";
		string1 = string1 + "<img src='images/icon_eps_15b.png' alt-text='[EPS] button' align='center' onclick='loadExternal(\""+mainFigure2+"eps\")'>";
		if (isFile(fig_linear_manual+"png")) {
			string1 = string1 + " : <font size='-1'>Linear-Manual</font> ";
			string1 = string1 + "<img src='images/icon_png_15b.png' alt-text='[PNG] button' align='center' onclick='loadImage(\""+key+"\",\""+fig_linear_manual+"png\",\"100\")'> ";
			string1 = string1 + "<img src='images/icon_eps_15b.png' alt-text='[EPS] button' align='center' onclick='loadExternal(\""+fig_linear_manual+"eps\")'>";
		}
		// GBrowse annotation file
		if (isFile(CGD_annotations_SNP)) {
			string1 = string1 + " : <font size='-1'>GBrowse-SNP</font> ";
			string1 = string1 + "<button onclick=\"loadExternal('"+CGD_annotations_SNP+"');\">GBrowse</button>";
		}
		if (isFile(url_base+"users/"+user+"/super.txt")) {  // Super-user privilidges.
			// Show CNV bias figure for WGseq and ddRADseq.
			if ((isFile(CNV_bias_WGseq_1)) || (isFile(CNV_bias_WGseq_2))) {
				string1 = string1 + " : <font size='-1'>CNV biases</font> ";
				if (isFile(CNV_bias_WGseq_1)) {
					string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+CNV_bias_WGseq_1+"\",\"100\")'>1</button>";
				} else if (isFile(CNV_bias_WGseq_2)) {
					string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+CNV_bias_WGseq_2+"\",\"100\")'>2</button>";
				}
			} else if ((isFile(CNV_bias_ddRADseq_1)) || (isFile(CNV_bias_ddRADseq_2)) || (isFile(CNV_bias_ddRADseq_3)) || (isFile(CNV_bias_ddRADseq_4))) {
				string1 = string1 + " : <font size='-1'>CNV biases</font> ";
				string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+CNV_bias_ddRADseq_1+"\",\"100\")'>1</button>";
				string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+CNV_bias_ddRADseq_2+"\",\"100\")'>2</button>";
				string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+CNV_bias_ddRADseq_3+"\",\"100\")'>3</button>";
				string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+CNV_bias_ddRADseq_4+"\",\"100\")'>T</button>";
			}
		}
		// Show allelic ratio plot version for ddRADseq.
		if ((isFile(fig_a_linear_SNPratio_ddRADseq)) || (isFile(fig_b_linear_SNPratio_ddRADseq)) || (isFile(fig_c_linear_SNPratio_ddRADseq)) || (isFile(fig_d_linear_SNPratio_ddRADseq))) {
			string1 = string1 + " : <font size='-1'>SNP ratios</font> ";
			if (isFile(fig_a_linear_SNPratio_ddRADseq)) {
				string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+fig_a_linear_SNPratio_ddRADseq+"\",\"100\")'>1</button>";
			}
			if (isFile(fig_b_linear_SNPratio_ddRADseq)) {
				string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+fig_b_linear_SNPratio_ddRADseq+"\",\"100\")'>2</button>";
			}
			if (isFile(fig_c_linear_SNPratio_ddRADseq)) {
				string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+fig_c_linear_SNPratio_ddRADseq+"\",\"100\")'>3</button>";
			}
			if (isFile(fig_d_linear_SNPratio_ddRADseq)) {
				string1 = string1 + "<button onclick='loadImage(\""+key+"\",\""+fig_d_linear_SNPratio_ddRADseq+"\",\"100\")'>4</button>";
			}
		}

		var string2 = "</td><td width='5%' align='right'><div onclick=\'closeProject(\""+user+"\",\""+project+"\",\""+key+"\",\""+color1+"\",\""+color2+"\",\""+parent+"\");' style='display:inline-block;'><b>[X]</b></div></td></tr>";
		var string3 = "<tr><td align='center' colspan='3'>";
		var string4 = "<div id='fig_"+key+"'><img src='"+mainFigure1+"png' width='100%'></div>";
		string4 = string4 + "</td></tr></table>"+"<hr></div>";

		visible_list.innerHTML += string1+string2+string3+string4;

		if (color1 != 'null') { // Hapmap analysis done.
			showColors(color1 ,'userProjectA_'+key,'a');
			showColors(color2 ,'userProjectB_'+key,'b');
			showColors('red' ,'userProjectHOM_'+key,'hom');
			showColors('grey','userProjectHET_'+key,'ab');
		} else { // No hapmap analysis done.
			if (parent == project) { // no LOH calculations done.
				showColors('grey','userProjectHET_'+key,'ab');
			} else { // LOH calculations have been done.
				showColors('red' ,'userProjectHOM_'+key,'hom');
				showColors('grey','userProjectHET_'+key,'ab');
			}
		}

		var show_button_element = document.getElementById("show_"+key);
		show_button_element.style.visibility='hidden';
		if (document.getElementById("delete_"+key)) {
			var delete_button_element = document.getElementById("delete_"+key);
			delete_button_element.style.visibility='hidden';
		}

		var projectsShown = localStorage.getItem("projectsShown");
		projectsShown = projectsShown+" "+user+":"+project+":"+key+":"+color1+":"+color2+":"+parent;
		localStorage.setItem("projectsShown", projectsShown);
	}
	function closeProject(user,project,key,color1,color2,parent) {
		var show_button_element = document.getElementById("show_"+key);
		show_button_element.style.visibility="visible";

		if(document.getElementById("delete_"+key)) {
			var delete_button_element = document.getElementById("delete_"+key);
			delete_button_element.style.visibility="visible";
		}

		var figure_element = document.getElementById("figure_"+key);
		figure_element.remove();

		var projectsShown = localStorage.getItem("projectsShown");
		projectsShown = projectsShown.replace(" "+user+":"+project+":"+key+":"+color1+":"+color2+":"+parent,"");
		localStorage.setItem("projectsShown", projectsShown);
	}
	</script>
<hr>
<div id="visible_list" name="visible_list"></div>

<script type="text/javascript">
<!---------------- Page consistency upon logout --------------!>
update_projectsShown_after_logout = function() {
	localStorage.removeItem('projectsShown');
}


<!---------------- Page consistency upon new project section --------------!>
update_projectsShown_after_new_project = function() {
	var projectsShown_entries = projectsShown.split(' ');
	var new_projectsShown = "";
	localStorage.setItem("projectsShown","");
	for (var i=0;i<projectsShown_entries.length; i++) {
		var currentProject = projectsShown_entries[i];
		if (currentProject != '') {
			var entry_parts    = currentProject.split(':');
			projID = parseInt(entry_parts[2])+1;
			projID.toString();
			new_projectsShown = new_projectsShown + entry_parts[0]+":"+entry_parts[1]+":"+projID+":"+entry_parts[3]+":"+entry_parts[4]+":"+entry_parts[5]+" ";
		}
	}
	new_projectsShown = new_projectsShown.slice(0, -1);
	localStorage.setItem("projectsShown",new_projectsShown);
}


<!---------------- Page consistency upon delete project section --------------!>
update_projectsShown_after_project_delete = function(deletedProjectKey) {
	// Clean up 'deletedProjectID' to projectID.   'p1_[number]_delete' => [number]
	var deletedProjectID = deletedProjectKey;
	deletedProjectID = deletedProjectID.replace("p1_","");
	deletedProjectID = deletedProjectID.replace("_delete","");
	deletedProjectID = Number(deletedProjectID);

	// adjust projectsShown string to reflect new positions of projects.
	var projectsShown_entries = projectsShown.split(' ');
	var new_projectsShown = "";
	localStorage.setItem("projectsShown","");
	for (var i=0;i<projectsShown_entries.length; i++) {
		var currentProject = projectsShown_entries[i];
		if (currentProject != '') {
			var entry_parts    = currentProject.split(':');
			projID = parseInt(entry_parts[2]);
			if (projID < deletedProjectID) {
				projID.toString();
				new_projectsShown = new_projectsShown + entry_parts[0]+":"+entry_parts[1]+":"+projID+":"+entry_parts[3]+":"+entry_parts[4]+":"+entry_parts[5]+" ";
			} else if (projID > deletedProjectID) {
				projID = projID-1;
				projID.toString();
				new_projectsShown = new_projectsShown + entry_parts[0]+":"+entry_parts[1]+":"+projID+":"+entry_parts[3]+":"+entry_parts[4]+":"+entry_parts[5]+" ";
			}
		}
	}
	new_projectsShown = new_projectsShown.slice(0, -1);
	localStorage.setItem("projectsShown",new_projectsShown);
}


<!---------------- Page consistency upon refresh section --------------!>
if(localStorage.getItem("tabInUse")){
	var tabInUse      = localStorage.getItem("tabInUse");
} else {
	var tabInUse      = "user";
}
if(localStorage.getItem("projectsShown")){
	var projectsShown = localStorage.getItem("projectsShown");
}
if (!projectsShown) {
	<?php
	$projectsShown = "";
	foreach ($systemProjectFolders as $key=>$project) {
		// Load colors for project.
		$colorFile        = $GLOBALS['directory']."users/default/projects/".$project."/colors.txt";
		if (file_exists($colorFile)) {
			$handle       = fopen($colorFile,'r');
			$colorString1 = trim(fgets($handle));
			$colorString2 = trim(fgets($handle));
			fclose($handle);
		} else {
			$colorString1 = 'null';
			$colorString2 = 'null';
		}
		$parentFile   = $GLOBALS['directory']."users/default/projects/".$project."/parent.txt";
		$handle       = fopen($parentFile,'r');
		$parentString = trim(fgets($handle));
		fclose($handle);
		$projectsShown = $projectsShown." default:".$project.":".($key+$userProjectCount).":null:null:".$parentString." ";
	}
	$projectsShown = trim($projectsShown);
	echo "var projectsShown = '".$projectsShown."';\n";
	?>
}
tabWindow_user();
<?php
if (isset($_SESSION['logged_on'])) {
	?>
	switch(tabInUse) {
		case "user":     tabWindow_user();     break
		case "project":  tabWindow_project();  break
		case "genome":   tabWindow_genome();   break
		case "hapmap":   tabWindow_hapmap();   break
		case "system":   tabWindow_system();   break
		case "about":    tabWindow_about();    break
		case "examples": tabWindow_examples(); break
	}
	<?php
} else {
	?>
	tabWindow_user();
	<?php
}
?>
// show projects.
if (projectsShown) {
	var projectsShown_entries = projectsShown.split(' ');
	localStorage.setItem("projectsShown","");
	for (var i=0;i<projectsShown_entries.length; i++) {
		var currentProject = projectsShown_entries[i];
		if (currentProject != '') {
//			// show string description in front of each projectList entry.
//			var entry_parts       = currentProject.split(':');
//			var projID            = parseInt(entry_parts[2])+1;
//			projID.toString();
//			var new_projectString = entry_parts[0]+":"+entry_parts[1]+":"+projID+":"+entry_parts[3]+":"+entry_parts[4]+":"+entry_parts[5];
//			var visible_list      = document.getElementById("visible_list");
//			visible_list.innerHTML += "["+new_projectString+"]<br>\n";

			// show figure for project.
			var entry_parts    = currentProject.split(':');
			openProject(entry_parts[0], entry_parts[1], entry_parts[2], entry_parts[3], entry_parts[4], entry_parts[5]);
		}
	}
}


<!---------------- Page reload when logged out --------------!>
function reload_page() {
	var go = false;
	// trigger reload if reload_page is called and user is found to be logged off.
	<?php
	if (!isset($_SESSION['logged_on'])) {
		echo "\t// User is not logged in.\n";
		//location.reload();
	} else {
		echo "// User is still logged in.\n";
	}
	?>
}
<?php
if (isset($_SESSION['logged_on'])) {
	echo "\n// Initiate recurrent call to reload_page function, which depends upon project status.\n";
	echo "var intervalID = window.setInterval(reload_page, 3000);\n";
}
?>
</script>
<!-- Manage interfaces for deleting objects.--!>
	<script type="text/javascript" src="js/project.delete_index.js"></script>
	<script type="text/javascript" src="js/genome.delete_index.js"></script>
	<script type="text/javascript" src="js/user.delete_index.js"></script>
	<script type="text/javascript" src="js/hapmap.delete_index.js"></script>
<!-- Needed for deletion functions.--!>
	<script type="text/javascript" src="js/jquery-1.7.2.js"></script>
	<script type="text/javascript" src="js/jquery.form.js"></script>

</body>
</html>
