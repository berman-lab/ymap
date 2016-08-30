function deleteHapmapConfirmation(user,hapmap,key){
	console.log("$ user      = '"+user+"'");
	console.log("$ project   = '"+hapmap+"'");
	console.log("$ key       = '"+key+"'");
	panel_iframe             = document.getElementById('panel_hapmap_iframe');
	dom_object               = panel_iframe.contentDocument.getElementById('h_delete_'+key);

	dom_object.innerHTML     = "<b><font color=\"red\">[Are you sure?]</font><button type='button' onclick='parent.deleteHapmap_yes(\""+user+"\",\""+hapmap+"\",\""+key+"\")'>Yes, delete.</button>";
	dom_object.innerHTML    += "<button type='button' onclick='parent.deleteHapmap_no(\""+user+"\",\""+hapmap+"\",\""+key+"\")'>No, cancel</button></b>";

	dom_button               = panel_iframe.contentDocument.getElementById('hapmap_delete_'+key);
	dom_button.style.display = 'none';
}

function deleteHapmap_yes(user,hapmap,key){
	$.ajax({
		url : 'hapmap.delete_server.php',
		type : 'post',
		data : {
			user: user,
			hapmap: hapmap
		},
		success : function(answer){
			if(answer == "COMPLETE"){
				// reload entire page - in order to ensure the update of the quota calculation
				window.top.location.reload(); 
			}
		}
	});
}


function deleteHapmap_no(user,hapmap,key){
	panel_iframe             = document.getElementById('panel_hapmap_iframe');
	dom_object               = panel_iframe.contentDocument.getElementById('h_delete_'+key)
	dom_object.innerHTML     = "";

	dom_button               = panel_iframe.contentDocument.getElementById('hapmap_delete_'+key);
	dom_button.style.display = 'inline';
}
