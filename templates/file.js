test = "ofijoij"

document.forms['data'].elements['input'].onchange = function(evt) {
    console.log("SODIJFSODIF");
    if(!window.FileReader) return; // Browser is not compatible

    var reader = new FileReader();

    reader.onload = function(evt) {
        if(evt.target.readyState != 2) return;
        if(evt.target.error) {
            alert('Error while reading file');
            return;
        }

        filecontent = evt.target.result;

        document.forms['data'].elements['fit_params'].value = evt.target.result;
		csvdata = evt.target.result;
    };

    reader.readAsText(evt.target.files[0]);
};
