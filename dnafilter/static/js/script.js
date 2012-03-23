$(function() {
  // focus on selected input method when radio changes
  $('input[name="input_method"]').change(function(e) {
    var target = $(e.target);
    $('*[name="' + target.val() + '"]').focus();
  });

  // update radio when input is added
  $('*[name="seqs"], *[name="seqdatafile"]').change(function(e) {
    var target = $(e.target);
    $('*[value="' + target.attr('name') + '"]').attr('checked', true);
  });

  // setup color-pickers
  $('.colorpicker').wheelColorPicker({
    dir: '/img/colorpicker'
  });
});
