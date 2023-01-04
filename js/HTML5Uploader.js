/*
 * Modern Uploader
 * https://github.com/filad/Modernuploader
 *
 * Copyright 2013, Adam Filkor
 * http://filkor.org
 *
 * Licensed under the MIT license:
 * http://www.opensource.org/licenses/MIT
 */

/*jslint nomen: true, regexp: true */
/*global define, window, URL, webkitURL, FileReader */

var allFiles = new Array;

(function (factory) {
	'use strict';
	if (typeof define === 'function' && define.amd) {	// Register as an anonymous AMD module:
		define(['jquery', 'tmpl', './jquery.fileupload-resize', './jquery.fileupload-validate'], factory);
	} else {						// Browser globals:
		factory(window.jQuery, window.Handlebars);
	}
}

(function ($, Handlebars) {
	'use strict';
	$.blueimp.fileupload.prototype._specialOptions.push( 'filesContainer', 'uploadTemplateId', 'downloadTemplateId' );

	// The UI version extends the file upload widget and adds complete user interface interaction:
	$.widget('filkor.html5Uploader', $.blueimp.fileupload, {
		options: {	autoUpload:		false,			// true = load on selection; false = load on upload button.
				maxChunkSize:		1024*512,		// break uploads into 0.5Mb fragments.
				uploadTemplateId:	'template-upload',	// The ID of the upload template.
				downloadTemplateId:	'template-download',	// The ID of the download template.
				filesContainer:		'#files',		// The file list container default name.
				dataType:		'json',			// The expected data type of the $.ajax upload request.
				getNumberOfFiles:	function(){		// Function to return the current number of files.   (used by maxNumberOfFiles validation.)
								return this.filesContainer.children().length;
							},
				add:			function (e, data) {	// CALLBACK : invoked when files are added to the fileupload widget.
								var $this     = $(this),
								that          = $this.data('filkor-html5Uploader'),
								options       = that.options,
								files         = data.files,
								addmore       = $("#add-more-button"),
								existingFiles = options.existingFiles || [];

								$.getJSON('php/uploader', {file: target_dir + data.files[0].name}, function (result) {
									var file = result.file;
									data.uploadedBytes = file && file.size;

									// Add file names to global scope variable 'allFiles'...   replacing comma characters with something not allowed in file names, to allow PHP to
									// split names by commas added by javascript.   PHP will replace the "\\" with "," later.
									allFiles.push(files[0].name.replace(",","\\"));

									data.process(function () {
										return $this.html5Uploader('process', data);
									}).always(function () {
									//	$('#info-wrapper-1').fadeIn(0);	// this section of uploader interface is handled in per-page defined scripts.
										$('#info-wrapper-2').fadeIn(0);
										$('#info-wrapper-3').fadeIn(0);
										$('#select-wrapper').fadeOut(0);
										addmore.removeClass('hidden');

										// 'data.context' represents a single 'li.file-item'; attached as an object to the 'data' attribute.
										data.context = that._renderUpload(files).data('data', data);
										options.filesContainer[ options.prependFiles ? 'prepend' : 'append' ](data.context);

										that._forceReflow(data.context);
										that._transition(data.context).done(
											function () {
												if ((that._trigger('added', e, data) !== false) && (options.autoUpload || data.autoUpload) && data.autoUpload !== false && !data.files.error) {
													data.submit();
												}
											}
										);
									});
								});
							},
				send:			function (e, data) {	// CALLBACK : invoked at the start of each file upload request.
								var that = $(this).data('filkor-html5Uploader');
								if (data.context && data.dataType && data.dataType.substr(0, 6) === 'iframe') {
									// Iframe Transport does not support progress events.   Lacking an indeterminate progress bar, we show the full animated bar.   This is a bit hacky.
									if (!$.support.transition) {
										data.context
											.find('.progress').hide()
											.find('.bar')
											.css('width','100%');
										data.context
											.find('.progress-wrap')
											.append('<img src="img/progress-static.gif" class="progress-static">')
											.css('height','7px');
									}
								}
								return that._trigger('sent', e, data);
							},
				done:			function (e, data) {	// CALLBACK : called upon completion of each file upload request.
								var that     = $(this).data('filkor-html5Uploader'),
								getFileList  = data.getFileList || that.options.getFileList,
								files        = getFileList(data),
								progressbar  = data.context.find('.progress .bar'),
								file         = files[0] || {error: 'Empty file upload result'};

								// could have used _transition, but it's buggy for some reason.
								data.context.find('.cancel-single').hide();

								if (!$.support.transition) {	//hide the 'static' progress animation
									data.context.find('.progress-static').hide();
									data.context.find('.progress').show();
								}
							},
				stop:			function (e, data) {	// CALLBACK : called upon conclusion of all uploads.
								// Hide timeer and speed counters.
								$('#info-wrapper-1').fadeOut(0);
								$('#info-wrapper-3').fadeOut(0);

								//==================================================================
								// Generate and submit a form to POST data to the concluding script.
								//------------------------------------------------------------------
								// Script to run at conclusion of data upload.
								var conclusion = document.createElement("form");
									conclusion.setAttribute("method","post");
									//conclusion.setAttribute("action","upload_processer.php");
									conclusion.setAttribute("action","scripts_seqModules/scripts_WGseq/project.single_WGseq.install_1.php");
								// add dataFormat to form, if defined.
								if ("project" in window) {
								//	console.log("[]"+);
									var input0 = document.createElement("input");
										input0.setAttribute("type","hidden");
										input0.setAttribute("name","dataFormat");
										input0.setAttribute("value",dataFormat);
										conclusion.appendChild(input0);
								}
								// add uploaded filename to form.
								var input1 = document.createElement("input");
									input1.setAttribute("type","hidden");
									input1.setAttribute("name","fileName");
									var outputFileNames = allFiles;
									input1.setAttribute("value",outputFileNames);
									conclusion.appendChild(input1);
								// add user to form.
								var input2 = document.createElement("input");
									input2.setAttribute("type","hidden");
									input2.setAttribute("name","user");
									input2.setAttribute("value",user);
									conclusion.appendChild(input2);
								// add genome or project to form.
								var input3 = document.createElement("input");
									input3.setAttribute("type","hidden");
									if ("genome" in window) {
										input3.setAttribute("name","genome");
										input3.setAttribute("value",genome);
									} else if ("project" in window) {
										input3.setAttribute("name","project");
										input3.setAttribute("value",project);
									}
									conclusion.appendChild(input3);
								// add key to form.
								var input4 = document.createElement("input");
									input4.setAttribute("type","hidden");
									input4.setAttribute("name","key");
									input4.setAttribute("value",key);
									conclusion.appendChild(input4);
								// append form to allow submit
								document.body.appendChild(conclusion);
								// Submit generated form : works in Safari, but not Firefox.
								conclusion.submit();
							},
				fail:			function (e, data) {	// CALLBACK : upon upload failure.
								data.context.each(function (index) {
								var file   = data.files[index];
									file.error = file.error || data.errorThrown || true;
									console.log(file.error);
								});
							},
				progress:		function (e, data) {	// CALLBACK : for updating the progress-bar during upload.
								if (data.context) {
									var progress = Math.floor(data.loaded / data.total * 100);
									data.context.find('.progress')
										.attr('aria-valuenow', progress)
										.find('.bar').css('width',progress + '%');
								}

								// Hides "Start Upload" buttons during progress.
								$('#info-wrapper-1').fadeOut(0);
							},
				progressall:		function (e, data) {	// CALLBACK : for updating the time and bitrate indications.
								var $this   = $(this),
								timeInfo    = $this.find('.time-info'),
								bitrateInfo = $this.find('.speed-info');
								timeInfo.find('span').html(	$this.data('filkor-html5Uploader')._renderTimeInfo(data)	);
								bitrateInfo.html(			$this.data('filkor-html5Uploader')._renderBitrateInfo(data)	);
							},
				processstart:		function () {
								//console.log('processstart..');
							},
				destroy:		function (e, data) {
								//destroy file.
								//By default when you click on the cancel buttonn you only abort the jqXHR, it doesn't delete the file.
								//(If you want to deletion  you can implement it here)
							},
				getFileList:		function (data) {	// CALLBACK : retrieve list of files from the server response.
								if (data.result && $.isArray(data.result.files)) {
									return data.result.files;
								}
								return [];
							}
			},
			_renderTemplate:	function (func, files) {
							if (!func) { return $(); }
							var result = func({
								files:   files,
								options: this.options
							});
							if (result instanceof $) { return result; }
							return $(this.options.templatesContainer).html(result).children();
						},
			_renderUpload:		function(files) {
							return this._renderTemplate(
								this.options.uploadTemplate,
								files
							);
						},
			_forceReflow:		function (node) {			// http://stackoverflow.com/questions/9016307/force-reflow-in-css-transitions-in-bootstrap
							return $.support.transition && node.length && node[0].offsetWidth;
						},
			_transition:		function (node) {
							var dfd = $.Deferred();
							if ($.support.transition) {
								node.on(
									$.support.transition.end,
									function (e) {
										// Make sure we don't respond to other transitions events in the container element, e.g. from button elements:
										if (e.target === node[0]) {
											node.unbind($.support.transition.end);
											dfd.resolveWith(node);
										}
									}
								);
							} else {
								dfd.resolveWith(node);
							}
							return dfd;
						},
			_formatBitrate:		function (bits) {
							if (typeof bits !== 'number') {	return ''; }
							if (bits >= 8589934592) {		return (bits / 1073741824 / 8).toFixed(2) + ' GB/s'; }
							if (bits >= 12388608) {			return (bits / 1048576 / 8).toFixed(1) + ' MB/s'; }
							if (bits >= 8192) {				return (bits / 1024 / 8).toFixed(0) + ' KB/s'; }
							if (bits < 0) {					return 0; }
							return (bits / 8).toFixed(2) + ' byte/s';
						},
			_formatTime:		function (seconds) {
							if (seconds < 0) { seconds = 0; }
							var date = new Date(seconds * 1000),
							days     = Math.floor(seconds / 86400);
							days     = days ? days + 'd ' : '';
							return days + ('0' + date.getUTCHours()).slice(-2) + ':' + ('0' + date.getUTCMinutes()).slice(-2) + ':' + ('0' + date.getUTCSeconds()).slice(-2);
						},
			_renderTimeInfo:	function (data) {
							return this._formatTime( (data.total - data.loaded) * 8 / data.bitrate );
						},
			_renderBitrateInfo:	function (data) {
							return this._formatBitrate(data.bitrate);
						},
			_initTemplates:		function () {		// Initialize & compile Handelbars template defined as id="template-upload".
							var options = this.options;
							options.templatesContainer = this.document[0].createElement(
								options.filesContainer.prop('nodeName')
							);
							if (Handlebars) {
								if (options.uploadTemplateId) {
									var source = $('#' + options.uploadTemplateId).html();
									options.uploadTemplate = Handlebars.compile(source);
								}
							}
						},
			_initFilesContainer:	function () {
							var options = this.options;
							if (options.filesContainer === undefined) {
								options.filesContainer = this.element.find('.files');
							} else if (!(options.filesContainer instanceof $)) {
								options.filesContainer = $(options.filesContainer);
							}
						},
	 	       _initHandlebarHelpers:	function () {
							//debug, usage {{debug}} or {{debug someValue}}
							Handlebars.registerHelper("debug", function (optionalValue) {
								console.log("Current Context");
								console.log("====================");
								console.log(this);

								if (optionalValue) {
									console.log("Value");
									console.log("====================");
									console.log(optionalValue);
								}
							});

							//format File size,
							Handlebars.registerHelper("formatFileSize", function (bytes) {
								if (typeof bytes !== 'number') {	return ''; }
								if (bytes >= 1073741824) {		return (bytes / 1073741824).toFixed(1) + ' GB'; }
								if (bytes >= 1048576) {			return (bytes / 1048576).toFixed(1) + ' MB'; }
													return (bytes / 1024).toFixed(0) + ' KB';
							});

							Handlebars.registerHelper("shortenName", function (name) {
								if (name.length > 45) { name = ' ' + name.substring(0, 45) + '...'; }
								return name;
							});
						},
			_startHandler:		function (e) {
							var data;
							e.preventDefault();
							$(".file-item").each(function(index, fileItem) {
								data = $(fileItem).data('data');
								if (data && data.submit && !data.jqXHR && !data.files.error && data.submit()) {
									//show pause btn for exmaple
								}
							});
						},
			_cancelHandler:		function (e) {		// CALLBACK : invoked when cancel button is pressed during upload.
							var template = $(e.currentTarget).closest('.file-item'),
							data         = template.data('data') || {},
							that         = this;

							// Reset user interface if file upload is interrrupted.
							$('#info-wrapper-1').fadeOut(0);
							$('#info-wrapper-2').fadeIn(0);
							$('#info-wrapper-3').fadeOut(0);
							$('#select-wrapper').fadeIn(0);

							template.slideUp('fast', function () {
								if (data.jqXHR) {
								data.jqXHR.abort();
									//we may also delete the file, even when it's partially uploaded
									//that._trigger('destroy', e, data);
								}
								template.remove();
							});
						},
			_initEventHandlers:	function () {
							var uploadBtn = $("#start-button");
							this._super();
							this._on(uploadBtn, {'click' : this._startHandler});
							this._on(this.options.filesContainer, {'click .cancel-single': this._cancelHandler} );
						},
			_initSpecialOptions:	function () {
							this._super();
							this._initFilesContainer();
							this._initTemplates();
						},
			_create:		function () {
							this._super();
							this._initHandlebarHelpers();
						}
	});
}));
