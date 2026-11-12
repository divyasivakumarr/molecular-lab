(function(root){
const radii={H:.31,C:.76,N:.71,O:.66,F:.57,P:1.07,S:1.05,Cl:1.02,Br:1.2,I:1.39};
function parse(text,name){let lines=text.replace(/\r/g,'').split('\n'),atoms=[],bonds=[],declared2D=false;
if(/V3000/.test(text))throw Error('Please export MOL V2000 or XYZ. V3000 is not supported.');
if(/V2000/.test(text)){let k=lines.findIndex(l=>l.includes('V2000')),n=+lines[k].slice(0,3),b=+lines[k].slice(3,6);declared2D=/2D/.test(lines[1]||'');if(!n||n>300)throw Error('Use a molecule with 1–300 atoms.');for(let i=0;i<n;i++){let l=lines[k+1+i]||'';atoms.push({e:l.slice(31,34).trim(),x:+l.slice(0,10),y:+l.slice(10,20),z:+l.slice(20,30)});}for(let i=0;i<b;i++){let l=lines[k+1+n+i]||'';bonds.push([+l.slice(0,3)-1,+l.slice(3,6)-1,+l.slice(6,9)]);}}
else if(/^\d+\s*$/.test(lines[0])){let n=+lines[0];if(!n||n>300)throw Error('Use a molecule with 1–300 atoms.');for(let i=0;i<n;i++){let a=(lines[i+2]||'').trim().split(/\s+/);atoms.push({e:a[0],x:Number(a[1]),y:Number(a[2]),z:Number(a[3])});}}
else throw Error('Upload a MOL V2000, SDF, or XYZ molecular file. Images and .moldraw project files cannot supply atom coordinates.');
if(atoms.some(a=>!radii[a.e]||![a.x,a.y,a.z].every(Number.isFinite)))throw Error('Invalid coordinates or unsupported element. Supported: H, C, N, O, F, P, S, Cl, Br, I.');
if(bonds.some(b=>b[0]<0||b[1]<0||b[0]>=atoms.length||b[1]>=atoms.length))throw Error('Invalid bond indices.');
let inferred=!/V2000/.test(text);if(inferred)for(let i=0;i<atoms.length;i++)for(let j=i+1;j<atoms.length;j++){let d=dist(atoms[i],atoms[j]);if(d>.35&&d<1.25*(radii[atoms[i].e]+radii[atoms[j].e]))bonds.push([i,j,1]);}
return {atoms,bonds,name,declared2D,flat:atoms.every(a=>Math.abs(a.z-atoms[0].z)<1e-6),inferred};}
function dist(a,b){return Math.hypot(a.x-b.x,a.y-b.y,a.z-b.z)}
function angle(a,b,c){let u=[a.x-b.x,a.y-b.y,a.z-b.z],v=[c.x-b.x,c.y-b.y,c.z-b.z],d=Math.hypot(...u)*Math.hypot(...v);return d?Math.acos(Math.max(-1,Math.min(1,u.reduce((s,x,i)=>s+x*v[i],0)/d)))*180/Math.PI:NaN;}
function measure(m){let lengths=m.bonds.map(([i,j])=>({label:`${m.atoms[i].e}${i+1}–${m.atoms[j].e}${j+1}`,value:dist(m.atoms[i],m.atoms[j])/10})),angles=[];m.atoms.forEach((a,c)=>{let ns=m.bonds.filter(b=>b[0]===c||b[1]===c).map(b=>b[0]===c?b[1]:b[0]);for(let i=0;i<ns.length;i++)for(let j=i+1;j<ns.length;j++)angles.push({label:`${m.atoms[ns[i]].e}${ns[i]+1}–${a.e}${c+1}–${m.atoms[ns[j]].e}${ns[j]+1}`,value:angle(m.atoms[ns[i]],a,m.atoms[ns[j]])});});return {lengths,angles};}
root.Chem={parse,dist,angle,measure};if(typeof module!=='undefined')module.exports=root.Chem;
})(globalThis);
