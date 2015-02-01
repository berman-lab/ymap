function setupAjaxForm(identifier){
	$(identifier).ajaxForm({
		beforeSubmit: function() {},
		success: showResponse
	});
}