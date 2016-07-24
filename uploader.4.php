<?php
	session_start();
?>
<!DOCTYPE html>
<!--[if lt IE 7]>      <html class="no-js lt-ie9 lt-ie8 lt-ie7"> <![endif]-->
<!--[if IE 7]>         <html class="no-js lt-ie9 lt-ie8"> <![endif]-->
<!--[if IE 8]>         <html class="no-js lt-ie9"> <![endif]-->
<!--[if gt IE 8]><!--> <html class="no-js"> <!--<![endif]-->
<head>
    <meta charset="UTF-8">
    <link rel="stylesheet" href="css/normalize.css">     <!-- Normalizes cross-browser display differences.  --!>
    <link rel="stylesheet" href="css/style.css">         <!-- Normalizes cross-browser styling.              --!>
    <link rel="stylesheet" href="css/defaultTheme.css">  <!-- --!>
	<style type="text/css">
		.tab {
			margin-left:    1cm;
		}
	</style>
</head>
<BODY class="tab">
<!---------- Main section of User interface ----------!>
	<form id="fileupload" class="HTML5Uploader" method="POST" action="uploader/" enctype="multipart/form-data">
		<div class="upload-wrapper">
			<div id="1-wrapper">
			<table><tr>
			<!---------- Add file to upload button. ----------!>
			<td valign="top">
				<div id="select-wrapper" class="info-wrapper">
					<div id="browsebutton" class="fileinput-button button gray" href="">
						<script type="text/javascript">
							console.log("uploader.4.php : display_string    = '"+display_string[0]+"'");
							console.log("uploader.4.php : currentDir        = '<?php echo getcwd(); ?>'");
							console.log("uploader.4.php : target_dir        = '"+target_dir+"'");
							console.log("uploader.4.php : conclusion_script = '"+conclusion_script+"'");
							document.write(display_string[0]);
						</script>
						<input type="file" id="fileinput" name="files[]" class="fileinput" single onchange="Show2()">
					</div>
				</div>
			</td>
			<!---------- Start upload button. ----------!>
			<td valign="top">
				<div id="info-wrapper-1" class="info-wrapper" style="display: none; font-size: 10px;">
					<button id="start-button" class="button greenish" type="submit">Upload</button>
				</div>
			</td>
			<!---------- Upload name and status. ----------!>
			<td valign="top">
				<div id="info-wrapper-2" class="info-wrapper" style="display: none; font-size: 10px;">
					<ul id="files" class="files">
				</div>
			</td>
			<!---------- Time and Speed of upload. ----------!>
			<td valign="top">
				<div id="info-wrapper-3" class="info-wrapper" style="display: none; font-size: 10px;">
					<div title="Remaining time" class="time-info"><span>00:00:00</span></div>
					<div title="Uploading speed" class="speed-info">0 KB/s</div>
				</div>
			</td>
			</tr></table>
			</div>

			<div id="2-wrapper" style="display:none">
			<table><tr>
			<!---------- Add file to upload button. ----------!>
			<td valign="top">
				<div id="select-wrapper" class="info-wrapper">
					<div id="browsebutton" class="fileinput-button button gray" href="">
						<script type="text/javascript">document.write(display_string[1]);</script>
						<input type="file" id="fileinput" name="files[]" class="fileinput" single onchange="Show3()">
					</div>
				</div>
			</td>
			<!---------- Start upload button. ----------!>
			<td valign="top">
				<div id="info-wrapper-1" class="info-wrapper" style="display: none; font-size: 10px;">
					<button id="start-button" class="button greenish" type="submit">Upload</button>
				</div>
			</td>
			<!---------- Upload name and status. ----------!>
			<td valign="top">
				<div id="info-wrapper-2" class="info-wrapper" style="display: none; font-size: 10px;">
					<ul id="files" class="files">
				</div>
			</td>
			<!---------- Time and Speed of upload. ----------!>
			<td valign="top">
				<div id="info-wrapper-3" class="info-wrapper" style="display: none; font-size: 10px;">
					<div title="Remaining time" class="time-info"><span>00:00:00</span></div>
					<div title="Uploading speed" class="speed-info">0 KB/s</div>
				</div>
			</td>
			</tr></table>
			</div>

			<div id="3-wrapper" style="display:none">
			<table><tr>
			<!---------- Add file to upload button. ----------!>
			<td valign="top">
				<div id="select-wrapper" class="info-wrapper">
					<div id="browsebutton" class="fileinput-button button gray" href="">
						<script type="text/javascript">document.write(display_string[2]);</script>
						<input type="file" id="fileinput" name="files[]" class="fileinput" single onchange="Show4()">
					</div>
				</div>
			</td>
			<!---------- Start upload button. ----------!>
			<td valign="top">
				<div id="info-wrapper-1" class="info-wrapper" style="display: none; font-size: 10px;">
					<button id="start-button" class="button greenish" type="submit">Upload</button>
				</div>
			</td>
			<!---------- Upload name and status. ----------!>
			<td valign="top">
				<div id="info-wrapper-2" class="info-wrapper" style="display: none; font-size: 10px;">
					<ul id="files" class="files">
				</div>
			</td>
			<!---------- Time and Speed of upload. ----------!>
			<td valign="top">
				<div id="info-wrapper-3" class="info-wrapper" style="display: none; font-size: 10px;">
					<div title="Remaining time" class="time-info"><span>00:00:00</span></div>
					<div title="Uploading speed" class="speed-info">0 KB/s</div>
				</div>
			</td>
			</tr></table>
			</div>

			<div id="4-wrapper" style="display:none">
			<table><tr>
			<!---------- Add file to upload button. ----------!>
			<td valign="top">
				<div id="select-wrapper" class="info-wrapper">
					<div id="browsebutton" class="fileinput-button button gray" href="">
						<script type="text/javascript">document.write(display_string[3]);</script>
						<input type="file" id="fileinput" name="files[]" class="fileinput" single onchange="Show5()">
					</div>
				</div>
			</td>
			<!---------- Start upload button. ----------!>
			<td valign="top">
				<div id="info-wrapper-1" class="info-wrapper" style="display: none; font-size: 10px;">
					<button id="start-button" class="button greenish" type="submit">Upload</button>
				</div>
			</td>
			<!---------- Upload name and status. ----------!>
			<td valign="top">
				<div id="info-wrapper-2" class="info-wrapper" style="display: none; font-size: 10px;">
					<ul id="files" class="files">
				</div>
			</td>
			<!---------- Time and Speed of upload. ----------!>
			<td valign="top">
				<div id="info-wrapper-3" class="info-wrapper" style="display: none; font-size: 10px;">
					<div title="Remaining time" class="time-info"><span>00:00:00</span></div>
					<div title="Uploading speed" class="speed-info">0 KB/s</div>
				</div>
			</td>
			</tr></table>
			</div>

			<!---------- Pass along the script to be run once all files are loaded. ----------!>
			<input type="hidden" id="hidden_field" name="target_dir" value="123">
			<script type="text/javascript">
				document.getElementById('hidden_field').value = target_dir;

			Show2=function() {
				// show second file select button.
				document.getElementById("2-wrapper").style.display = 'inline';

				// hide upload button until all files are selected.
				document.getElementById("info-wrapper-1").style.display = 'none';
			}
			Show3=function() {
				// hide second file select button.
				document.getElementById("2-wrapper").style.display = 'none';
				// show third file select button.
				document.getElementById("3-wrapper").style.display = 'inline';
			}
			Show4=function() {
				// hide third file select button.
				document.getElementById("3-wrapper").style.display = 'none';
				// show fourth file select button.
				document.getElementById("4-wrapper").style.display = 'inline';
			}
			Show5=function() {
				// hide fourth file select button.
				document.getElementById("4-wrapper").style.display = 'none';

				// show upload button once all files are selected.
				document.getElementById("info-wrapper-1").style.display = 'inline';
			}
			</script>
		</div>
	</form>

<!---------- Script section ----------!>
	<script id="template-upload" type="text/x-handlebars-template">
		{{#each files}}
			<li class="file-item">
				<div class="first">
					<span class="top">
						<span class="filename">{{shortenName name}}</span>
					</span>
					{{#if error}}
					<div class="error">
						<span>Error: </span>
						<span>{{error}}</span>
					</div>
					{{/if}}
				</div>
				<div class="second">
					<span class="filesize">{{formatFileSize size}}</span>
					<div class="progress-wrap">
						<div class="progress">
							<div class="bar" style="width:0%"></div>
						</div>
					</div>
					<div class="clear"></div>
				</div>
			</li>
		{{/each}}
	</script>

	<script src="http://ajax.googleapis.com/ajax/libs/jquery/1.9.1/jquery.min.js"></script>
	<script>window.jQuery || document.write('<script src="js/jquery-1.9.1.min.js"><\/script>')</script>

	<!-- handlebars -->
	<script src="js/handlebars.min.js"></script>

	<!-- HTML5Uploader -->
	<script src="js/jquery.ui.widget.js"></script>
	<script src="js/bootstrap-transition.js"></script>
	<script src="js/jquery.iframe-transport.js"></script>
	<script src="js/jquery.fileupload.js"></script>
	<script src="js/jquery.fileupload-process.js"></script>
	<script src="js/jquery.fileupload-validate.js"></script>
	<script src="js/HTML5Uploader.js"></script>
	<script src="js/main.js"></script>

</body>
</html>
