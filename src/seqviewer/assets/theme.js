(function () {
  try {
    var stored = localStorage.getItem('__STORAGE_KEY__');
    if (stored === 'dark') {
      document.documentElement.setAttribute('data-theme', 'dark');
    }
  } catch (e) {}
})();
