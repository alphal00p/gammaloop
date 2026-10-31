"""Independent exact-method selection; no wrappers around content, reads, or assertions."""
import hashlib, json, os, sys, unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, '/build/source/crates/linnet-py/tests')
import test_svg_browser as product
import symbolica.core as core
from symbolica.community import hepkit, feynkit_native
from playwright.sync_api import BrowserType, Browser
from playwright._impl import _transport

assert __debug__
assert os.environ['LINNET_TEST_BROWSER_ENGINES'] == 'chromium,webkit'
assert hashlib.sha256(Path(product.__file__).read_bytes()).hexdigest() == 'ee700e0629fad62358e1f0a778353757a9a305bdbc131f2ae2675b7b74d82c6c'
assert core.__file__.startswith('/work/venv/')
assert hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest() == 'd4cbe4e0fe79ebe128d8729b0d7081cb6eb2de1aee575d568fc533301098d957'
for cls in (hepkit.SquaredAmplitude, feynkit_native.SquaredAmplitude):
    assert all(hasattr(cls,n) for n in ('paged','_display_','_repr_mimebundle_'))
draws, engines, pages, subtests = [], [], [], []
original_draw, original_launch = product.pinned_drawing, BrowserType.launch
original_page, original_env = Browser.new_page, _transport.get_driver_env
private = ('SYMBOLICA_LICENSE','SYMBOLICA_LICENSE_PATH','SYMBOLICA_LICENSE_SIGNED','SYMBOLICA_LICENSE_FILE')
driver_calls = []

def driver_env():
    environment = original_env()
    for key in private: environment.pop(key,None)
    assert all(k not in environment for k in private)
    driver_calls.append(True)
    return environment

def launch(owner,*args,**kwargs):
    environment = dict(kwargs.get('env') or os.environ)
    for key in private: environment.pop(key,None)
    kwargs['env'] = environment
    browser = original_launch(owner,*args,**kwargs)
    engines.append({'engine':owner.name,'version':browser.version,'default_executable':owner.executable_path,'requested_executable':kwargs.get('executable_path'),'private_variables_absent':all(k not in environment for k in private)})
    return browser

def draw(x,y):
    result = original_draw(x,y)
    draws.append([x,y])
    return result

def new_page(owner,*args,**kwargs):
    page = original_page(owner,*args,**kwargs)
    # Only blank-page DPR is read. No extra action follows product set_content.
    pages.append({'dpr':page.evaluate('devicePixelRatio'),'viewport':page.viewport_size})
    return page

class Result(unittest.TextTestResult):
    def addSubTest(self,test,subtest,err):
        subtests.append({'engine':subtest.params.get('browser'),'passed':err is None})
        super().addSubTest(test,subtest,err)

suite = unittest.TestSuite([product.SvgBrowserTests('test_readable_scale_and_fitted_drawings_are_not_cropped')])
with patch.object(_transport,'get_driver_env',driver_env), patch.object(BrowserType,'launch',launch), patch.object(Browser,'new_page',new_page), patch.object(product,'pinned_drawing',draw):
    result = unittest.TextTestRunner(verbosity=2,resultclass=Result).run(suite)
complete = result.testsRun == 1 and result.wasSuccessful() and not result.skipped and draws == [[0,2],[22,2],[0,42],[64,2]] and subtests == [{'engine':'chromium','passed':True},{'engine':'webkit','passed':True}] and len(pages) == 2
record = {'complete':complete,'tests_run':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),'skips':len(result.skipped),'drawings':draws,'engines':engines,'pages':pages,'subtests':subtests,'driver_environment_calls':len(driver_calls),'source_sha256':hashlib.sha256(Path(product.__file__).read_bytes()).hexdigest(),'installed_core_sha256':hashlib.sha256(Path(core.__file__).read_bytes()).hexdigest(),'installed_core_path':core.__file__,'native_display_methods_present':True,'wrappers_on_page_content_locator_reads_or_assertions':False}
Path('/work/product.json').write_text(json.dumps(record,indent=2)+'\n')
assert complete, record
