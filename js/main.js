$(function () {
	'use strict';

	// Set up our file upload script on server.
	// defines some variables in dom element "fileupload" (form defined in uploader.1.php).
	$("#fileupload").html5Uploader({
		url: 'php/uploader/index.php',
		maxFileSize: 1024*1024*20 // 20MB
	});

	// Load existing files and set user direcory.
	$.ajax({
		url: $('#fileupload').html5Uploader('option', 'url'),
		dataType: 'json',
		context: $('#fileupload')[0]
	}).always(function () {

	}).done(function (result) {
		// Do with the result what you want here; this will return the already uploaded files (even if they're partial).
		//console.log(result.files);
	});
});
