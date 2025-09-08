#rm -rf GL_QQX_AAA; ../../target/dev-optim/gammaloop -s GL_QQX_AAA -o -t generation run generate_qqx_aaa.yaml
../../target/dev-optim/gammaloop -s GL_QQX_AAA -n -t inspect inspect --process-id 0 --name qqx_aaa_subtracted --discrete-dim 0 -m -p 0.1 0.2 0.3
rm -rf ./integration_workspace; ../../target/dev-optim/gammaloop -s GL_QQX_AAA -n -t integrate integrate --process-id 0 --name qqx_aaa_subtracted --workspace-path ./integration_workspace --result-path ./integration_workspace/integration_results.txt -c 8
