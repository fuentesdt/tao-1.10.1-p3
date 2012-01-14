To run minsurf example:
ccafe-single --ccafe-rc CCAFERC



To run minsurf with gui:
1. copy runOneProcWGUI script from ccaffeine tree to local directory
2. replace --ccafe-rc $pkgdatadir/CcaffeineRC with --ccafe-rc CCAFERC_gui
3. if desired, add --buildFile CCAFERC_gui.bld to $gui --builderPort $UIPORT
3. ./runOneProcWGUI
4. Use the 'configure' button on the TaoSolver component to change the
solver parameters
5. push 'go' button on TaoDriver component






